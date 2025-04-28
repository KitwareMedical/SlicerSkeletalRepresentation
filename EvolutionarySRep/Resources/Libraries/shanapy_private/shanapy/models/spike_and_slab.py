import numpy as np
import os
import scipy
def perturb(d):
    perturbation = 1e-4
    vInv = 2 * perturbation * np.eye(d)
    return vInv
def sigmoid(x):
    return 1 / (1 + np.exp(-x))
GRADIENT_MIN = 1e-6
STEP_MIN = 1/(2 ** 8)

class SpikeNSlab():
    """
    Training spike and slab prior for linking structure
    
    """
    def __init__(self, train_X=None, train_y=None, dim_per_link=14, p_all_hat=0.1, output_file=None, redo_selection=True):
        """
        Args:
            pos_data (numpy.ndarray): d-by-n design matrix, where n is the # of positive cases
            neg_data (numpy.ndarray): d-by-n design matrix, where n is the # of negative cases
            dim_per_link (int): the # of dimension of every link in every case
        """

        self.pos_data = train_X[:, train_y]
        self.neg_data = train_X[:, ~train_y]

        self.d, self.num_pos = self.pos_data.shape
        self.num_links = self.d // dim_per_link
        self.dim_per_link = dim_per_link

        self.redo_selection = redo_selection
        ## prior of a link being selected
        self.p_all_hat = p_all_hat
        self.output_file = output_file
        if os.path.exists(output_file) and redo_selection:
            os.system("rm " + output_file)

    def link_feature_ids(self, link_id):
        feature_id_start = link_id * self.dim_per_link
        feature_id_end = (link_id + 1) * self.dim_per_link - 1
        return feature_id_start, feature_id_end
    def update_zs(self, X, zs, t_hat, vs, B, p_all):
        """
        Update zs to maximize the likelihood: P(z | X, theta)
        Refer paper: Generalized Spike-and-Slab Priors for Bayesian Group Feature
Selection Using Expectation Propagation Appendix A.

        Args:
            X: design matrix
            zs (vector): current zs, indicators for every link, that need to be updated
            t_hat (vector): estimated noisy labels (classification) / values (regression)
            vs (scalar): 'slab' prior
            B (matrix):
            p_all (scalar): p_alpha_l in the paper
        Returns:
            L_llh_all (scalar): sum of log likelihood
            zs (vector): updated indicators
        """
        L_llh_all = 0.0
        unit_t_hat = t_hat / np.linalg.norm(t_hat)
        llh_for_links = np.zeros((1,self.num_links))
        A_minus_g = vs * np.eye(self.d)
        ## exclude unselected link features
        for unselected_link_id in np.where(zs==0)[0]:
            id_start, id_end = self.link_feature_ids(unselected_link_id)
            A_minus_g[id_start:id_end+1] = 0

        for i in range(self.num_links):
            L_total = .0
            if zs[i] == 1:
                ## exclude current link features
                curr_id_start, curr_id_end = self.link_feature_ids(i)
                ## minus this group (i.e., features related to this link)
                A_minus_g[curr_id_start:curr_id_end+1] = 0
                C_minus_g = B + X.T @ A_minus_g @ X
                ## Features of this group (i.e., link)
                X_g = X[curr_id_start:curr_id_end+1, :]
                C_inv = np.linalg.inv(C_minus_g + perturb(C_minus_g.shape[0]))
                M = X_g @ C_inv @ X_g.T

                eig_val, eig_vec = np.linalg.eig(M)
                L = .0
                for j in range(len(eig_val)):
                    s_j = eig_val[j]
                    e_j = eig_vec[:, j]
                    temp = unit_t_hat @ C_inv @ X_g.T
                    q_j = temp @ e_j
                    L = L + 0.5 * (q_j**2/(s_j + 1/vs) - np.log(1 + vs * s_j))
                L_total = np.exp(L)

                L_total = 1 if np.isnan(L_total) or np.isinf(L_total) else L_total

            else:
                L_total = 1.0
            llh = p_all * L_total / (p_all * L_total + (1-p_all))
            L_llh_all += np.log(llh)
            zs[i] = 1 if np.random.rand() < llh or L_total > 1000000000000.0 else 0
        return L_llh_all, zs
    def update_others(self, X, zs, vs, ws, selected_feat_ids, unselected_feat_ids):
        """
        Update other parameters given updated zs
        """
        ## update scalar p_all
        ## TODO: tune p_all_hat
        num_ones = len(np.where(zs==1)[0])
        num_zeros = self.num_links - num_ones
        p_all = np.random.beta(self.num_links * self.p_all_hat + num_ones, num_zeros+self.num_links * (1-self.p_all_hat))
        ## update vs
        ## TODO: tune parameters here
        vs = 1 / np.random.gamma(num_ones/2 + 1.5, 1)
        A = vs * np.eye(self.d)
        A_remain = A[selected_feat_ids, :]#np.delete(A, unselected_feat_ids, 0)
        A_remain = A_remain[:, selected_feat_ids]#np.delete(A_remain, unselected_feat_ids, 1)
        ws_nonzero = ws[selected_feat_ids]
        return p_all, vs, ws_nonzero, A_remain
    def optimize_ws(self, X_nonzero, ws_nonzero, ws, a_diag, y, selected_feat_ids, unselected_feat_ids):
        ### For debug
        # zs, ws, y, X, t_hat, vs, B, load_vars = self.load_mat_file('/playpen/workspace/newuoa/Correspondence_Oct7/shapeStat_upSide/Correspondence/test_param.mat')

        # X_nonzero, ws_nonzero, a_diag = load_vars['X_nz'], load_vars['ws_nz'].squeeze(),  np.diag(load_vars['A'])
        # selected_feat_ids = []
        # for selected_link_id in np.where(zs==1)[0]:
        #     feat_start_id, feat_end_id = self.link_feature_ids(selected_link_id)
        #     feat_ids = [k for k in range(feat_start_id, feat_end_id+1)]
        #     selected_feat_ids += feat_ids
        # unselected_feat_ids = []
        # for unselected_link_id in np.where(zs==0)[0]:
        #     feat_start_id, feat_end_id = self.link_feature_ids(selected_link_id)
        #     feat_ids = [k for k in range(feat_start_id, feat_end_id+1)]
        #     unselected_feat_ids += feat_ids
        ####
        errs = []
        curr_err, y_hat_temp = self._compute_error(np.dot(ws_nonzero, X_nonzero), y)
        max_iter_irls = 50

        for iw in range(max_iter_irls):
            errs.append(curr_err)
            e = y - y_hat_temp
            g = np.dot(X_nonzero, e) - (a_diag * ws_nonzero)
            if np.all(abs(g) < GRADIENT_MIN):
                ## current parameters are optimum
                break

            ## Hessian matrix
            sig = sigmoid(np.dot(ws_nonzero, X_nonzero))
            b_diag = sig * (1 - sig)
            B = np.diag(b_diag)
            H = X_nonzero @ B @ X_nonzero.T + np.diag(a_diag)
            upperH = np.linalg.cholesky(H + perturb(H.shape[0])).T
            step = 1
            while step > STEP_MIN:
                temp = ws_nonzero - np.dot(np.linalg.inv(step * upperH), np.dot(np.linalg.inv(upperH.T), g))
                data_error, y_hat_temp = self._compute_error(np.dot(temp, X_nonzero), y)
                regularizer = np.dot(a_diag.T, temp ** 2) / 2
                temp_total_error = data_error;
                if temp_total_error >= np.min(errs):
                    step /= 2
                else:
                    ws_nonzero = temp

                    new_total_error = temp_total_error
                    step = 0
        inv_upperH = np.linalg.inv(upperH)
        sigma = inv_upperH @ inv_upperH.T
        sig = sigmoid(np.dot(ws_nonzero, X_nonzero))
        b_diag = sig * (1-sig)
        B = np.diag(b_diag) + perturb(len(b_diag))
        t_hat = np.dot(ws_nonzero, X_nonzero) + np.dot(np.linalg.inv(B), y-y_hat_temp)

        sqrt_S = scipy.linalg.sqrtm(sigma)
        b_lap = sqrt_S / np.sqrt(2)
        ## sample ws_nonzero from laplacian
        def laprnd(mu, b_lap, size):
            u = np.random.rand(*size) - 0.5
            return mu - np.dot(b_lap, np.sign(u)) * np.log(1 - 2 * np.abs(u))
        temp_ws_nonzero = laprnd(ws_nonzero, b_lap, (len(ws_nonzero), ))

        # temp_p_ws = np.dot(np.exp(np.dot(-np.abs(temp_ws_nonzero - ws_nonzero), np.linalg.inv(b_lap))), np.linalg.inv(2 * b_lap))
        # p_ws[selected_feat_ids] = temp_p_ws
        # p_ws[unselected_feat_ids] = 0
        ws_nonzero = temp_ws_nonzero

        ws = np.zeros(self.d)
        ws[selected_feat_ids] = ws_nonzero

        return ws_nonzero, ws, t_hat, np.min(errs)

    def _compute_error(self, pred, target):
        """
        binary cross entropy loss:
        sum(y*log(y) + (1-y) * log(1-y))
        """
        y = sigmoid(pred)
        total_loss = .0
        eps = 1e-7
        for i, t in enumerate(target):
            if t == 1:
                curr_loss = np.log(eps) if y[i] < eps else np.log(y[i])
                total_loss += curr_loss
            else:
                curr_loss = np.log(eps) if (1-y[i]) < eps else np.log(1-y[i])
                total_loss += curr_loss

        if np.any(np.iscomplex(y)) or np.any(np.iscomplex(total_loss)):
            print("Complex loss was found")
        return -total_loss, y
    # def load_mat_file(self, file_name='/playpen/workspace/Simulate_Shapes/shanapy/models/test_python_sns.mat'):
    #     from scipy import io
    #     new_vars = io.loadmat(file_name)
    #     self.d = len(new_vars['ws'].squeeze())
    #     return new_vars['zs'].squeeze(), new_vars['ws'].squeeze(), new_vars['y'].squeeze(), new_vars['X'], new_vars['t_hat'].squeeze(), new_vars['vs'].squeeze(), new_vars['B'], new_vars
    def load_zs(self):

        with open(self.output_file, 'rb') as f:
            train_err = np.load(f)
            zs = np.load(f)
            cum_zs = np.load(f)
        return train_err, cum_zs, zs
    def fit(self, pred=False, test_X=None, test_y=None):
        if not self.redo_selection: return self.load_zs()
        _, num_pos = self.pos_data.shape
        d, num_neg = self.neg_data.shape
        X = np.concatenate([self.pos_data, self.neg_data], axis=1) # d-by-N, where N = #pos + #neg
        #### 1. Initialize variables
        ## initialize all links to be selected (p = 1)
        zs = np.random.binomial(1, 1, (self.num_links,))
        y = np.concatenate([np.ones((1, num_pos)), np.zeros((1, num_neg))], axis=1).squeeze()
        ws = np.random.normal(size=(d, ))
        t_hat = y
        vs = 0.06
        B = np.eye(X.shape[1])
        p_all = 0.1
        ws_nonzero = None
        selected_link_ids = None
        ce_errs = []
        min_err = 2000
        best_zs = []

        cum_zs = np.zeros_like(zs)
        burn_in = 1000
        # zs, ws, y, X, t_hat, vs, B, _ = self.load_mat_file()
        mcmc_samples = burn_in * 2
        ## MCMC sampling process
        for t in range(1, mcmc_samples + burn_in + 1):
            #################
            ##  1. update zs
            #################
            L_llh, zs = self.update_zs(X, zs, t_hat, vs, B, p_all)

            if np.all([i==0 for i in zs]):
                ## no links were selected
                print('No links were selected.')
                continue
            selected_feat_ids = []
            for selected_link_id in np.where(zs==1)[0]:
                feat_start_id, feat_end_id = self.link_feature_ids(selected_link_id)
                feat_ids = [k for k in range(feat_start_id, feat_end_id+1)]
                selected_feat_ids += feat_ids
            unselected_feat_ids = []
            for unselected_link_id in np.where(zs==0)[0]:
                feat_start_id, feat_end_id = self.link_feature_ids(unselected_link_id)
                feat_ids = [k for k in range(feat_start_id, feat_end_id+1)]
                unselected_feat_ids += feat_ids

            #################
            ##  2. update other hyperparameters
            #################
            p_all, vs, ws_nonzero, A = self.update_others(X, zs, vs, ws, selected_feat_ids, unselected_feat_ids)
            ### Now estimate/refine Weights all ws
            X_nonzero = X[selected_feat_ids]

            #################
            ##  3. update weights of classifier
            #################
            ws_nonzero, ws, t_hat, err = self.optimize_ws(X_nonzero, ws_nonzero, ws, np.diag(A), y, selected_feat_ids, unselected_feat_ids)
            ### TODO: record the process
            train_err, y_hat_train = self._compute_error(np.dot(ws_nonzero, X_nonzero), y)
            ce_errs.append(train_err)
            if np.any(np.iscomplex(train_err)):
                print("Complex error was found")
            if t > burn_in:
                cum_zs += zs

            if err < min_err:
                best_zs = zs
                min_err = err
            test_err = np.inf
            if pred:
                selected_test_feats = test_X[selected_feat_ids, :]
                test_err, y_hat_test = self._compute_error(np.dot(ws_nonzero, selected_test_feats), test_y)
                print('Iteration %d, Training loss: %f, Test loss: %f, selected %d links' % (t, train_err, test_err, len(np.where(zs==1)[0])))
            else:
                print('iteration %d, err: %f, selected %d links' % (t, train_err, len(np.where(zs==1)[0])))

        last_selection = np.where(best_zs==1)[0]
        with open(self.output_file, 'wb') as f:
            np.save(f, np.array(ce_errs))
            np.save(f, last_selection)
            np.save(f, np.array(cum_zs))
            # f.write("%s" % train_err)
            # if pred:
            #     f.write(",%s" % test_err)
            # for selected_link_id in np.where(zs==1)[0]:
            #     f.write(",%s" % selected_link_id)
            # f.write("\n")

        return ce_errs, cum_zs, last_selection
    def fit_and_pred(self, test_X, test_y):
        self.fit(True, test_X, test_y)
# model = SpikeNSlab()
# model.fit()
