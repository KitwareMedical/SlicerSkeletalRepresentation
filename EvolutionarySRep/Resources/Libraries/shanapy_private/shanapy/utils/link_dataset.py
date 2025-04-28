import os
from torch.utils.data import Dataset
import numpy as np
class LinkDataSet(Dataset):
    def __init__(self, data_dir=None):
        """
        Data are n-by-d matrices
        """
        ## n-by-d
        self.pos, self.neg = self._read(data_dir)

        self.data = np.concatenate([self.pos, self.neg], axis=0)
        self.y = np.concatenate([np.ones(self.pos.shape[0]), np.zeros(self.neg.shape[0])])
    def __len__(self):
        if self.data is None: return 0
        return self.data.shape[0]
    def __getitem__(self, idx):
        if self.data is None: return None
        return self.data[idx, :]
    def _read(self, data_dir=None):
        pos_dir = os.path.join(data_dir, 'pos')
        neg_dir = os.path.join(data_dir, 'neg')

        pos_data = self._parse_links(pos_dir)
        neg_data = self._parse_links(neg_dir)

        return pos_data, neg_data
    def _parse_links(self, data_dir=None):
        ret = [] # n-by-d matrix
        for case_file in os.listdir(data_dir):
            #print('Read linking ' + case_file + '...')
            file_name = os.path.join(data_dir, case_file)
            links = [] # k-by-14 matrix
            with open(file_name, 'r') as f:
                for str_spoke_tuple in f.readlines():
                    spoke_tuple = str_spoke_tuple.split(';')
                    link_feature = [] # 1-by-14 vector
                    link_feature += [float(x) for x in spoke_tuple[0].split(',')]
                    link_feature += [float(x) for x in spoke_tuple[1].split(',')]
                    links.append(link_feature)

            np_links = np.array(links)
            ret.append(np_links.flatten())
        return np.array(ret)
# loader = DataLoader('/playpen/workspace/my_paper/linking/data/links/')
