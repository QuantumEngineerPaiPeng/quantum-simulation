"""
This part of code is the DQN brain, which is a brain of the agent.
All decisions are made in here.
Using Tensorflow to build the neural network.

View more on my tutorial page: https://morvanzhou.github.io/tutorials/

Using:
Tensorflow: 1.0
gym: 0.8.0
"""

import numpy as np
import pandas as pd
import tensorflow as tf
from copy import deepcopy


# Deep Q Network off-policy
class DeepQNetwork:
    def __init__(
            self,
            n_actions,
            n_features,
            maxTime,
            learning_rate=0.00025,
            reward_decay=0.9,
            e_greedy=0.9,
            replace_target_iter=300,
            memory_size=500,
            batch_size=32,
            e_greedy_increment=None,
            output_graph=False,
            end_action=5,
            node_each_layer=[10]
    ):
        self.n_actions = n_actions
        self.n_features = n_features
        self.lr = learning_rate
        self.gamma = reward_decay
        self.epsilon_max = e_greedy
        self.replace_target_iter = replace_target_iter
        self.memory_size = memory_size
        self.batch_size = batch_size
        self.epsilon_increment = e_greedy_increment
        self.epsilon = 0 if e_greedy_increment is not None else self.epsilon_max

        # consist of [target_net, evaluate_net]
        self.node_each_layer = node_each_layer
        self.node_each_layer.insert(0, n_features)

        if output_graph:
            # $ tensorboard --logdir=logs
            # tf.train.SummaryWriter soon be deprecated, use following
            tf.summary.FileWriter("logs/", self.sess.graph)

        self.end_action = end_action
        self.maxTime = maxTime

        self.reset()

    def reset(self):
        # total learning step
        self.learn_step_counter = 0
        self.memory_counter = 0
        self.elite_memory_counter = 0

        # initialize zero memory [s, a, r, s_]
        self.memory = np.zeros((self.memory_size, self.n_features * 2 + 2))
       
        self.elite_memory = np.zeros((self.memory_size, self.n_features * 2 + 2))

        tf.reset_default_graph()  # added by Chao
        self._build_net()
        t_params = tf.get_collection('target_net_params')
        e_params = tf.get_collection('eval_net_params')
        self.replace_target_op = [tf.assign(t, e) for t, e in zip(t_params, e_params)]

        self.sess = tf.Session()

        self.sess.run(tf.global_variables_initializer())
        self.cost_his = []

        self.epsilon = 0

    def _build_net(self):
        # ------------------ build evaluate_net ------------------
        self.s = tf.placeholder(tf.float32, [None, self.n_features], name='s')  # input
        self.q_target = tf.placeholder(tf.float32, [None, self.n_actions], name='Q_target')  # for calculating loss
        with tf.variable_scope('eval_net'):
            # c_names(collections_names) are the collections to store variables
            c_names, w_initializer, b_initializer = \
                ['eval_net_params', tf.GraphKeys.GLOBAL_VARIABLES], \
                tf.random_normal_initializer(0., 0.3), tf.constant_initializer(0.1)  # config of layers

            l_list = [self.s for l_index in range(len(self.node_each_layer))]
            for layer_i in range(len(self.node_each_layer)-1):
                with tf.variable_scope('l'+str(layer_i)):
                    w1 = tf.get_variable('w'+str(layer_i), [self.node_each_layer[layer_i], \
                        self.node_each_layer[layer_i+1] ], initializer=w_initializer, collections=c_names)
                    b1 = tf.get_variable('b'+str(layer_i), [1, self.node_each_layer[layer_i+1]], \
                        initializer=b_initializer, collections=c_names)
                    l_list[layer_i+1] = tf.nn.relu(tf.matmul(l_list[layer_i], w1) + b1)

            with tf.variable_scope('l_out'):
                w_out = tf.get_variable('w_out', [self.node_each_layer[-1], self.n_actions], \
                    initializer=w_initializer, collections=c_names)
                b_out = tf.get_variable('b_out', [1, self.n_actions], initializer=b_initializer, collections=c_names)
                self.q_eval = tf.matmul(l_list[-1], w_out) + b_out
            '''
            # 2 layers
            # first layer. collections is used later when assign to target net
            with tf.variable_scope('l1'):
                w1 = tf.get_variable('w1', [self.n_features, n_l1], initializer=w_initializer, collections=c_names)
                b1 = tf.get_variable('b1', [1, n_l1], initializer=b_initializer, collections=c_names)
                l1 = tf.nn.sigmoid(tf.matmul(self.s, w1) + b1)

            # second layer. collections is used later when assign to target net
            with tf.variable_scope('l2'):
                w2 = tf.get_variable('w2', [n_l1, self.n_actions], initializer=w_initializer, collections=c_names)
                b2 = tf.get_variable('b2', [1, self.n_actions], initializer=b_initializer, collections=c_names)
                self.q_eval = tf.matmul(l1, w2) + b2
            '''

        with tf.variable_scope('loss'):
            self.loss = tf.reduce_mean(tf.squared_difference(self.q_target, self.q_eval))
        with tf.variable_scope('train'):
#             self._train_op = tf.train.RMSPropOptimizer(self.lr,momentum=0.95,epsilon=0.01).minimize(self.loss)
            self._train_op = tf.train.AdamOptimizer(self.lr).minimize(self.loss)

        # ------------------ build target_net ------------------
        self.s_ = tf.placeholder(tf.float32, [None, self.n_features], name='s_')    # input
        with tf.variable_scope('target_net'):
            # c_names(collections_names) are the collections to store variables
            c_names = ['target_net_params', tf.GraphKeys.GLOBAL_VARIABLES]

            l_list = [self.s_ for l_index in range(len(self.node_each_layer))]
            for layer_i in range(len(self.node_each_layer)-1):
                with tf.variable_scope('l'+str(layer_i)):
                    w1 = tf.get_variable('w'+str(layer_i), [self.node_each_layer[layer_i], \
                        self.node_each_layer[layer_i+1] ], initializer=w_initializer, collections=c_names)
                    b1 = tf.get_variable('b'+str(layer_i), [1, self.node_each_layer[layer_i+1]], \
                        initializer=b_initializer, collections=c_names)
                    l_list[layer_i+1] = tf.nn.relu(tf.matmul(l_list[layer_i], w1) + b1)

            with tf.variable_scope('l_out'):
                w_out = tf.get_variable('w_out', [self.node_each_layer[-1], self.n_actions], \
                    initializer=w_initializer, collections=c_names)
                b_out = tf.get_variable('b_out', [1, self.n_actions], initializer=b_initializer, collections=c_names)
                self.q_next = tf.matmul(l_list[-1], w_out) + b_out

    def store_transition(self, s, a, r, s_):
        if not hasattr(self, 'memory_counter'):
            self.memory_counter = 0

        transition = np.hstack((s, [a, r], s_))

        # replace the old memory with new memory
        index = self.memory_counter % self.memory_size
        self.memory[index, :] = transition

        self.memory_counter += 1
        
    def elite_store_transition(self, s, a, r, s_):
        if not hasattr(self, 'elite_memory_counter'):
            self.elite_memory_counter = 0

        transition = np.hstack((s, [a, r], s_))

        # replace the old memory with new memory
        index = self.elite_memory_counter % self.memory_size
        self.elite_memory[index, :] = transition

        self.elite_memory_counter += 1

    def choose_action(self, observation, epsilon=None):
        # to have batch dimension when feed into tf placeholder
        observation = observation[np.newaxis, :]
        epsilon = self.epsilon if epsilon==None else epsilon
        if np.random.uniform() < self.epsilon:
            # forward feed the observation and get q value for every actions
            actions_value = self.sess.run(self.q_eval, feed_dict={self.s: observation})
            action = np.argmax(actions_value)
        else:
            action = np.random.randint(0, self.n_actions)
        return action

    def learn(self, alpha=None, ifelite=False, endgame=False):
        if alpha == 0:
            return
        # check to replace target parameters
        if self.learn_step_counter % self.replace_target_iter == 0:
            self.sess.run(self.replace_target_op)
            #print('\ntarget_params_replaced\n')
        if not ifelite:
            # sample batch memory from all memory
            if self.memory_counter > self.memory_size:
                sample_index = np.random.choice(self.memory_size, size=self.batch_size)
            else:
                sample_index = np.random.choice(self.memory_counter, size=self.batch_size)
            batch_memory = self.memory[sample_index, :]
        else:
            # sample batch memory from elite memory
            if self.elite_memory_counter > self.memory_size:
                sample_index = np.random.choice(self.memory_size, size=self.batch_size)
            else:
                sample_index = np.random.choice(self.elite_memory_counter, size=self.batch_size)
            batch_memory = self.elite_memory[sample_index, :]

        q_next, q_eval = self.sess.run(
            [self.q_next, self.q_eval],
            feed_dict={
                self.s_: batch_memory[:, -self.n_features:],  # fixed params
                self.s: batch_memory[:, :self.n_features],  # newest params
            })

        # change q_target w.r.t q_eval's action
        #q_target = q_eval.copy()
        q_target = deepcopy(q_eval)

        batch_index = np.arange(self.batch_size, dtype=np.int32)
        eval_act_index = batch_memory[:, self.n_features].astype(int)
        state_n_index = batch_memory[:, -1].astype(int)
        reward = batch_memory[:, self.n_features + 1]

        # enggame should be treated differently in order not to blow up, chao yin
        for batch_index_i in batch_index:
            if eval_act_index[batch_index_i]  == self.end_action or (1 in batch_memory[batch_index_i][-6::1] ):
                q_target[batch_index_i, eval_act_index[batch_index_i]] = reward[batch_index_i]
            else:
                q_next_max = self.gamma * np.max(q_next, axis=1) 
                q_target[batch_index_i, eval_act_index[batch_index_i]] = reward[batch_index_i] + q_next_max[batch_index_i]
        """
        For example in this batch I have 2 samples and 3 actions:
        q_eval =
        [[1, 2, 3],
         [4, 5, 6]]

        q_target = q_eval =
        [[1, 2, 3],
         [4, 5, 6]]

        Then change q_target with the real q_target value w.r.t the q_eval's action.
        For example in:
            sample 0, I took action 0, and the max q_target value is -1;
            sample 1, I took action 2, and the max q_target value is -2:
        q_target =
        [[-1, 2, 3],
         [4, 5, -2]]

        So the (q_target - q_eval) becomes:
        [[(-1)-(1), 0, 0],
         [0, 0, (-2)-(6)]]

        We then backpropagate this error w.r.t the corresponding action to network,
        leave other action as error=0 cause we didn't choose it.
        """

        # train eval network
        if alpha != None:
            temp = self.lr
            self.lr = alpha
        _, self.cost = self.sess.run([self._train_op, self.loss],
                                     feed_dict={self.s: batch_memory[:, :self.n_features],
                                                self.q_target: q_target})
        if alpha != None:
            self.lr = temp

        self.cost_his.append(self.cost)

        # increasing epsilon
        self.epsilon = self.epsilon + self.epsilon_increment if self.epsilon < self.epsilon_max else self.epsilon_max
        self.learn_step_counter += 1

    def plot_cost(self):
        import matplotlib.pyplot as plt
        plt.plot(np.arange(len(self.cost_his)), self.cost_his)
        plt.ylabel('Cost')
        plt.xlabel('training steps')
        plt.show()



