
import numpy as np


class Markov_Chain:
    def __init__(self, state, initial, transition, step):
        # state_name_list is 1 by n array of int
        # initial_prob_dist is 1 by n array of floats
        # transition_matrix is n by n array of floats
        # step is the number of simulated state in a markov chain
        self.initial = initial
        self.transition = transition
        self.step = step
        self.state = state

    def simulator(self):
        # array index starts from ZERO!!

        # how many states are there
        # P is m by n array, len(P)=m
        num_of_state = len(self.initial)

        # initial an empty array to fill state index, 1 by step
        # range(3)=0,1,2
        final = [0 for i in range(self.step)]

        # array of cumulative initial probabilities, 1 by n+1
        c_initial = [0 for i in range(num_of_state + 1)]

        for i in range(num_of_state):
            # P[0][0:1]=P[0][0];
            c_initial[i + 1] = sum(self.initial[0:i + 1])

        # simulate the first state
        first = np.random.uniform(0, 1)
        for i in range(num_of_state):
            if (first >= c_initial[i]) & (first <= c_initial[i + 1]):
                temp = i

        # now using transition matrix to simulate the rest of the chain

        # matrix of cumulative transition matrix
        c_transition = [[0 for i in range(num_of_state + 1)]
                        for j in range(num_of_state)]

        for i in range(num_of_state):
            for j in range(num_of_state):
                c_transition[i][j + 1] = sum(self.transition[i][0:j + 1])

        # loop to update each state in the chain
        for i in range(self.step):
            final[i] = temp
            prob = c_transition[temp]
            first = np.random.uniform(0, 1)
            for j in range(num_of_state):
                if (first >= prob[j]) & (first <= prob[j + 1]):
                    temp = j

        # translate "final"(the state index) to state
        chain = [0 for i in range(self.step)]
        for i in range(self.step):
            chain[i] = self.state[final[i]]

        return chain

    def invariant_dist(self):
        # P=[1,2,3], R=asmatrix(P) is a 1 by 3 matrix,len(P)=3, len(R)=1
        # Q=[[1],[2],[3]], S=asmatrix(Q) is a 3 by 1 matrix, len(Q)=3, len(S)=3
        # matrix index can only use [i,j] instead of [i][j] for array

        # transfer array into matrix
        m_initial = np.asmatrix(self.initial)
        m_transition = np.asmatrix(self.transition)

        # matrix multiplication notice:
        # P*Q^2=(P*Q)^2 != P*Q*Q

        # track the number of m_transition in product_old
        mtply = 1
        product_old = m_initial * m_transition
        product_new = product_old * m_transition

        while (np.any(abs(product_new - product_old) > 0.000000000001)):
            product_old = product_new
            mtply = mtply + 1
            product_new = product_new * m_transition

        print(mtply)
        return product_old
