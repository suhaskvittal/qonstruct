
import networkx as nx

import ldpc

def title(name):
    print('######## %s ########' % name)

title('GIRTH TEST')

gr = ldpc.make_regular_tanner_graph(3, 4, 5)
print(nx.girth(gr))

title('PARITY CHECK MATRIX')

H = ldpc.tanner_graph_to_parity_check_matrix(gr)
print(H.astype(int))
