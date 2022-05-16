#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import random
import argparse
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import torch
import torch_geometric

from cell_type_annotation_model import DNNModel, SpatialModelTrainer


random.seed(0)
np.random.seed(0)
torch.manual_seed(0)
torch.cuda.manual_seed(0)
torch.cuda.manual_seed_all(0)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True


config = {
    'data': {
        'data_dir': 'dataset/MERFISH/',
        'save_dir': 'model/MERFISH/',
        'dataset': 'MERFISH',
    },
    'preprocess': {
        'filter_mt': True,
        'cell_min_counts': 300,
        'gene_min_cells': 10,
        'cell_max_counts_percent': 98.0,
        'drop_rate': 0,
    },
    'transfer': {
        'dnn_model': 'dnn_model/checkpoint_MERFISH_s.t7',
        'gpu': '0',
        'batch_size': 4096,
    },
    'train': {
        'pca_dim': 200,  # for BGI only
        'k_graph': 30,
        'edge_weight': True,
        'kd_T': 1,
        'feat_dim': 64,
        'w_dae': 1.0,
        'w_gae': 1.0,
        'w_cls': 10.0,
        'epochs': 200,
    }
}


def spatial_classification_tool(config, data_name):
    ''' Spatial classification workflow.

    # Arguments
        config (Config): Configuration parameters.
        data_name (str): Data name.
    '''
    ######################################
    #         Part 1: Load data          #
    ######################################
    
    # Set path and load data.
    print('\n==> Loading data...')
    dataset = config['data']['dataset'] 
    data_dir, save_dir = config['data']['data_dir'], config['data']['save_dir'] 
    print(f'  Data name: {data_name} ({dataset})')
    print(f'  Data path: {data_dir}')
    print(f'  Save path: {save_dir}')
    adata = sc.read_h5ad(os.path.join(data_dir, f'{data_name}.h5ad'))

    # Initalize save path.
    model_name = f'modelGDAE-{data_name}'
    save_dir = os.path.join(save_dir, model_name)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)


    ######################################
    #         Part 2: Preprocess         #
    ######################################

    print('\n==> Preprocessing...')
    strings = [f'{k}={v}' for k, v in config['preprocess'].items()]
    print('  Parameters(%s)' % (', '.join(strings)))

    # Preprocess data.
    if dataset == 'BGI':
        params = config['preprocess']
        if params['filter_mt']:
            adata.var['mt'] = adata.var_names.str.startswith(('Mt-', 'mt-'))
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
            adata = adata[adata.obs['pct_counts_mt'] < 10].copy()
        if params['cell_min_counts'] > 0:
            sc.pp.filter_cells(adata, min_counts=params['cell_min_counts'])
        if params['gene_min_cells'] > 0:
            sc.pp.filter_genes(adata, min_cells=params['gene_min_cells'])
        if params['cell_max_counts_percent'] < 100:
            max_counts = np.percentile(adata.X.sum(1), params['cell_max_counts_percent'])
            sc.pp.filter_cells(adata, max_counts=max_counts)
    if type(adata.X) != np.ndarray:
        adata_X_sparse_backup = adata.X.copy()
        adata.X = adata.X.toarray()
    print('  %s: %d cells × %d genes.' % (data_name, adata.shape[0], adata.shape[1]))

    # Please be aware: 
    #     DNN model takes the origin gene expression matrix through its own normalization as input.
    #     Other normalization (e.g. scanpy) can be added after DNN model inference is completed.

    # Add noise manually.
    if dataset != 'BGI':
        drop_factor = (np.random.random(adata.shape) > config['preprocess']['drop_rate']) * 1.
        adata.X = adata.X * drop_factor


    ######################################
    #  Part 3: Transfer from sc-dataset  #
    ######################################

    print('\n==> Transfering from sc-dataset...')
    strings = [f'{k}={v}' for k, v in config['transfer'].items()]
    print('  Parameters(%s)' % (', '.join(strings)))

    # Set device.
    os.environ['CUDA_VISIBLE_DEVICES'] = config['transfer']['gpu']
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load DNN model trained by sc-dataset.
    checkpoint = torch.load(config['transfer']['dnn_model'])
    dnn_model = checkpoint['model'].to(device)
    dnn_model.eval()

    # Initialize DNN input.
    marker_genes = checkpoint['marker_genes']
    gene_indices = adata.var_names.get_indexer(marker_genes)
    adata_X = np.pad(adata.X, ((0,0),(0,1)))[:, gene_indices].copy()
    norm_factor = np.linalg.norm(adata_X, axis=1, keepdims=True)
    norm_factor[norm_factor == 0] = 1
    dnn_inputs = torch.Tensor(adata_X / norm_factor).split(config['transfer']['batch_size'])

    # Inference with DNN model.
    dnn_predictions = []
    with torch.no_grad():
        for batch_idx, inputs in enumerate(dnn_inputs):
            inputs = inputs.to(device)
            outputs = dnn_model(inputs)
            dnn_predictions.append(outputs.detach().cpu().numpy())
    label_names = checkpoint['label_names']
    adata.obsm['pseudo_label'] = np.concatenate(dnn_predictions)
    adata.obs['pseudo_class'] = pd.Categorical([label_names[i] for i in adata.obsm['pseudo_label'].argmax(1)])
    adata.uns['pseudo_classes'] = label_names

    # Compute accuracy (only for MERFISH).
    if dataset == 'MERFISH':
        indices = np.where(~adata.obs['subclass'].isin(['L4/5 IT', 'L6 IT Car3', 'PVM', 'other']))[0]
        adjusted_pr = adata.obs['pseudo_class'][indices].to_numpy()
        adjusted_gt = adata.obs['subclass'][indices].replace(['Micro'], ['Macrophage']).to_numpy()
        acc = (adjusted_pr == adjusted_gt).sum() / len(indices) * 100.0
        print('  %s Acc (transfer only): %.2f%%' % (data_name, acc))


    ######################################
    #      Part 4: Train GDAE model      #
    ######################################

    print('\n==> Model training...')
    strings = [f'{k}={v}' for k, v in config['train'].items()]
    print('  Parameters(%s)' % (', '.join(strings)))

    # Normalize gene expression.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata_X = (adata.X - adata.X.mean(0)) / (adata.X.std(0) + 1e-10)

    # Construct spatial graph.
    gene_mat = torch.Tensor(adata_X)
    if dataset == 'BGI':  # PCA
        u, s, v = torch.pca_lowrank(gene_mat, config['train']['pca_dim'])
        gene_mat = torch.matmul(gene_mat, v)
    cell_coo = torch.Tensor(adata.obsm['spatial'])
    data = torch_geometric.data.Data(x=gene_mat, pos=cell_coo)
    data = torch_geometric.transforms.KNNGraph(k=config['train']['k_graph'], loop=True)(data)
    data.y = torch.Tensor(adata.obsm['pseudo_label'])

    # Make distances as edge weights.
    if config['train']['edge_weight']:
        data = torch_geometric.transforms.Distance()(data)
        data.edge_weight = 1 - data.edge_attr[:,0]
    else:
        data.edge_weight = torch.ones(data.edge_index.size(1))

    # Train self-supervision model.
    input_dim = data.num_features
    num_classes = len(adata.uns['pseudo_classes'])
    trainer = SpatialModelTrainer(input_dim, num_classes, device, config['train'])
    trainer.train(data, config['train'])
    trainer.save_checkpoint(os.path.join(save_dir, f'{model_name}.t7'))

    # Inference.
    print('\n==> Inferencing...')
    predictions = trainer.valid(data)
    celltype_pred = pd.Categorical([adata.uns['pseudo_classes'][i] for i in predictions.argmax(1)])
    if dataset == 'MERFISH':
        indices = np.where(~adata.obs['subclass'].isin(['L4/5 IT', 'L6 IT Car3', 'PVM', 'other']))[0]
        adjusted_pr = celltype_pred[indices].to_numpy()
        adjusted_gt = adata.obs['subclass'][indices].replace(['Micro'], ['Macrophage']).to_numpy()
        acc = (adjusted_pr == adjusted_gt).sum() / len(indices) * 100.0
        print('  %s Acc (transfer+GDAE): %.2f%%' % (data_name, acc))

    # Save results.
    result = pd.DataFrame({'cell': adata.obs_names.tolist(), 'celltype_pred': celltype_pred})
    result.to_csv(os.path.join(save_dir, f'{model_name}.csv'), index=False)
    adata.obsm['celltype_prob'] = predictions
    adata.obs['celltype_pred'] = pd.Categorical(celltype_pred)
    if 'adata_X_sparse_backup' in locals():
        adata.X = adata_X_sparse_backup
    adata.write(os.path.join(save_dir, f'{model_name}.h5ad'))

    # Save visualization.
    spot_size = (30 if dataset == 'BGI' else 20)
    if dataset == 'BGI':
        pseudo_top100 = adata.obs['pseudo_class'].to_numpy()
        other_classes = list(pd.value_counts(adata.obs['pseudo_class'])[100:].index)
        pseudo_top100[adata.obs['pseudo_class'].isin(other_classes)] = '_Others'
        adata.obs['pseudo_class'] = pd.Categorical(pseudo_top100)
    sc.pl.spatial(adata, img_key=None, color=['pseudo_class'], spot_size=spot_size, show=False)
    plt.savefig(os.path.join(save_dir, f'pseudo-{data_name}.pdf'), bbox_inches='tight', dpi=150)
    sc.pl.spatial(adata, img_key=None, color=['celltype_pred'], spot_size=spot_size, show=False)
    plt.savefig(os.path.join(save_dir, f'{model_name}.pdf'), bbox_inches='tight', dpi=150)
    print('  Predictions is saved in', os.path.join(save_dir, f'{model_name}.csv/pdf'))


if __name__ == '__main__':
    data_list = ['mouse1_sample1', 'mouse1_sample2', 'mouse1_sample3', 'mouse1_sample4', 'mouse1_sample5', 'mouse1_sample6',
                 'mouse2_sample1', 'mouse2_sample2', 'mouse2_sample3', 'mouse2_sample4', 'mouse2_sample5', 'mouse2_sample6']
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_name', choices=data_list)
    args = parser.parse_args()
    spatial_classification_tool(config, args.data_name)
