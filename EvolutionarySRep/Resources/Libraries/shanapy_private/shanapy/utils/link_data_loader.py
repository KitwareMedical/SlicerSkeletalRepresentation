import os
import numpy as np
class DataLoader():
    def __init__(self, data_dir=None):
        """
        Data are n-by-d matrices
        """
        self.pos, self.neg = self._read(data_dir)
    def _read(self, data_dir=None):
        pos_dir = os.path.join(data_dir, 'pos_link')
        neg_dir = os.path.join(data_dir, 'neg_link')

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
