# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 15:28:43 2020

Author: Rounak Meyur
Description: Test program to implement Deep Reinforcement Learning for Frozen Lake 
game.
"""

import sys
import numpy as np
import gym
import random

env = gym.make("Taxi-v3")
env.render()

action_size = env.action_space.n
state_size = env.observation_space.n

qtable = np.zeros((state_size, action_size))

## Hyperparameter settings
total_episodes = 50000        # Total episodes
learning_rate = 0.8           # Learning rate
max_steps = 99                # Max steps per episode
gamma = 0.95                  # Discounting rate

# Exploration parameters
epsilon = 1.0                 # Exploration rate
max_epsilon = 1.0             # Exploration probability at start
min_epsilon = 0.01            # Minimum exploration probability 
decay_rate = 0.005            # Exponential decay rate for exploration prob

## Q-learning procedure
rewards = []
for episode in range(total_episodes):
    state = env.reset()
    step = 0
    done = False
    total_rewards = 0
    
    for step in range(max_steps):
        exp_exp_tradeoff = random.uniform(0,1)
        if exp_exp_tradeoff > epsilon: action = np.argmax(qtable[state,:])
        else: action = env.action_space.sample()
        
        new_state, reward, done, info = env.step(action)
        qtable[state, action] = qtable[state, action] + \
            learning_rate * (reward + gamma * np.max(qtable[new_state, :]) - \
                             qtable[state, action])
        total_rewards += reward
        state = new_state
        if done == True:
            break
    
    # Update epsilon for new episode
    epsilon = min_epsilon + (max_epsilon - min_epsilon)*np.exp(-decay_rate*episode) 
    rewards.append(total_rewards)

print ("Score over time: " +  str(sum(rewards)/total_episodes))
print(qtable)

sys.exit(0)
#%% Play the game
env.reset()
rewards = []
for episode in range(100):
    state = env.reset()
    step = 0
    done = False
    total_rewards = 0
    print("****************************************************")
    print("EPISODE ", episode+1)

    for step in range(max_steps):
        
        # Take the action (index) that have the maximum expected future reward given that state
        action = np.argmax(qtable[state,:])
        
        new_state, reward, done, info = env.step(action)
        
        total_rewards += reward
        
        if done:
            rewards.append(total_rewards)
            
            # We print the number of step it took.
            print("Score", total_rewards)
            break
        state = new_state
env.close()