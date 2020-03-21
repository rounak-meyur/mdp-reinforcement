# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 15:28:43 2020

Author: Rounak Meyur
Description: Test program to implement Deep Reinforcement Learning for Inverted
Pendulum Problem.
"""

import gym

env = gym.make('SpaceInvaders-v0')
env.reset()

for _ in range(1000):
    env.render()
    env.step(env.action_space.sample())
    
env.close()