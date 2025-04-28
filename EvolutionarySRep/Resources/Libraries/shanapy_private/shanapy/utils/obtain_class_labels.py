import pickle
import os
def generate_class_labels(path, groups=['pos_refined', 'neg_refined'], output_path='../../data/class_labels.pkl'):
    pos_ids = []
    neg_ids = []
    pos_path = os.path.join(path, groups[0])
    neg_path = os.path.join(path, groups[1])
    for case_id in os.listdir(pos_path):
        pos_ids.append(case_id)
    for case_id in os.listdir(neg_path):
        neg_ids.append(case_id)

    with open(output_path, 'wb') as f:
        pickle.dump([pos_ids, neg_ids],f)
def load_class_labels(path='../../data/', source_path='/playpen/ra_job/jiyao/hippocampus/'):
    """
    Return pos class ids and neg class ids
    """
    file_name = os.path.join(path, 'class_labels.pkl')
    if not os.path.exists(file_name):
        generate_class_labels(source_path)
    with open(file_name, 'rb') as f:
        pos_ids, neg_ids = pickle.load(f)
    return pos_ids, neg_ids

if __name__ == '__main__':
    pos, neg = load_class_labels(source_path='/playpen/ra_job/jiyao/hippocampus/')
    print("Positive: ")
    print(pos)
    print("Negative: ")
    print(neg)