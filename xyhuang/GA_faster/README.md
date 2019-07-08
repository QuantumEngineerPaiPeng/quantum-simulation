# Genetic Algorithm

arXiv:1712.06567v3 [cs.NE] 20 Apr 2018\
https://github.com/paraschopra/deepneuroevolution.git\

with mutation only\
Possible way out: Crossover, fully random agents...\

NN: CartPoleAI() start with state dimension; end with action dimension\

Env initialization: run_agents()\

mutation power = 0.8 (able to converge) 0.005(unstable, but it is actually appropriate to use because the range of  parameter in NN is about this-----> 0.8 is too large )\


number of agents in one generation =2000 or could be even larger say 20000\
parents select = 20 only consider_top = 10, should be larger if number of agents is increased\

I only tried 100 generations, but it actually converges quite fast(mutation power = 0.8). I am not sure whether (mutation power = 0.005) could get some other results by going beyond 100 generation.



