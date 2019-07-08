from gym.envs.registration import register

register(
    id='pp-v0',
    entry_point='gym_pp.envs:PpEnv',
)