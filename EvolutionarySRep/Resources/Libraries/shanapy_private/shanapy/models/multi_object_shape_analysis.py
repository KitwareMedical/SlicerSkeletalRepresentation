### This file collect more spokes for the multi-object shape analysis
### Deprecated
# import vtk
# import os
# import logging
# import numpy as np
# # from spoke import Spoke
# # from s2_geometry import *

# def recon_sreps(dirs):
#     """
#     According to new dirs (from AJIVE) reconstruct s-reps

#     Return newly reconstructed s-reps
#     """
#     pass

# def diff_dirs(result, ground_truth):
#     """
#     Compare difference between result and ground truth

#     Input result: n x d
#     Input ground_truth: n x d

#     Return histogram of cos_sim

#     """
#     n, d = result.shape
#     ret = []
#     num_spokes = d // 3
#     for i in range(n):
#         dirs_mat = np.reshape(result[i, :], (num_spokes, 3)) # k x 3
#         gt_dirs_mat = np.reshape(ground_truth[i, :], (num_spokes, 3)) # k x 3

#         for j in range(num_spokes):
#             cos_sim = np.dot(dirs_mat[j, :], gt_dirs_mat[j, :])
#             ret.append(cos_sim)
#     return np.array(ret)
# def dim_reduct(A):
#     """
#     Dimension reduction according to screeplots
#     """
#     U, S, V = np.linalg.svd(A)
#     eigvals = S**2 / np.sum(S**2)  # NOTE (@amoeba): These are not PCA eigenvalues. 
#                                    # This question is about SVD.

#     from matplotlib import pyplot as plt
#     import matplotlib

#     fig = plt.figure(figsize=(8,5))
#     sing_vals = np.arange(len(eigvals)) + 1
#     plt.plot(sing_vals, eigvals, 'ro-', linewidth=2)
#     plt.title('Scree Plot')
#     plt.xlabel('Principal Component')
#     plt.ylabel('Eigenvalue')
#     #I don't like the default legend so I typically make mine like below, e.g.
#     #with smaller fonts and a bit transparent so I do not cover up data, and make
#     #it moveable by the viewer in case upper-right is a bad place for it 
#     leg = plt.legend(['Eigenvalues from SVD'], loc='best', borderpad=0.3, 
#                      shadow=False, prop=matplotlib.font_manager.FontProperties(size='small'),
#                      markerscale=0.4)
#     leg.get_frame().set_alpha(0.4)
#     leg.draggable(state=True)
#     plt.show()