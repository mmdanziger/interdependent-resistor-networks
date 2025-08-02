# interdependent-resistor-networks
Code for "Interdependent resistor networks with process-based dependency"

This code contains the simulation code including voltage calculations on the disordered lattice and various plotter/helper functions in python. There is also a python implementation of the basic logic in [./rrun.py](./rrun.py).

When this code was written in 2014, there were no adequate Python matrix calculations so the hard work was done in C++. Maybe the situation is different today.

If you use this code, please cite:

```bibtex
@article{Danziger_2015,
doi = {10.1088/1367-2630/17/4/043046},
url = {https://dx.doi.org/10.1088/1367-2630/17/4/043046},
year = {2015},
month = {apr},
publisher = {IOP Publishing},
volume = {17},
number = {4},
pages = {043046},
author = {Danziger, Michael M and Bashan, Amir and Havlin, Shlomo},
title = {Interdependent resistor networks with process-based dependency},
journal = {New Journal of Physics},
abstract = {Studies of resilience of interdependent networks have focused on structural dependencies between pairs of nodes across networks but have not included the effects of dynamic processes taking place on the networks. Here we study the effect of dynamic process-based dependencies on a system of interdependent resistor networks. We describe a new class of dependency in which a node’s functionality is determined by whether or not it is actually carrying current and not just by its structural connectivity to a spanning component. This criterion determines its functionality within its own network as well as its ability to provide support-but not electrical current-to nodes in another network. We present the effects of this new type of dependency on the critical properties of σ and , the overall conductivity of the system and the fraction of nodes which carry current, respectively. Because the conductance of current has direct physical effects (e.g. heat, magnetic induction), the development of a theory of process-based dependency can lead to innovative technology. As an example, we describe how the theory presented here could be used to develop a new kind of highly sensitive thermal or gas sensor.}
}
```

This repo also includes code for the exponentially distributed link lengths research, which was concomitant in my PhD. Specifically the exp_cc.bin binary that is built by [main_topo.cpp](./idresnet/main_topo.cpp). If you use that part please cite:


```bibtex
@article{Danziger_2016,
doi = {10.1209/0295-5075/115/36002},
url = {https://dx.doi.org/10.1209/0295-5075/115/36002},
year = {2016},
month = {sep},
publisher = {EDP Sciences, IOP Publishing and Società Italiana di Fisica},
volume = {115},
number = {3},
pages = {36002},
author = {Danziger, Michael M. and Shekhtman, Louis M. and Berezin, Yehiel and Havlin, Shlomo},
title = {The effect of spatiality on multiplex networks},
journal = {Europhysics Letters},
abstract = {Many multiplex networks are embedded in space, with links more likely to exist between nearby nodes than distant nodes. For example, interdependent infrastructure networks can be represented as multiplex networks, where each layer has links among nearby nodes. Here, we model the effect of spatiality on the robustness of a multiplex network embedded in 2-dimensional space, where links in each layer are of variable but constrained length. Based on empirical measurements of real-world networks, we adopt exponentially distributed link lengths with characteristic length ζ. By changing ζ, we modulate the strength of the spatial embedding. When ζ → ∞, all link lengths are equally likely, and the spatiality does not affect the topology. However, when  only short links are allowed, and the topology is overwhelmingly determined by the spatial embedding. We find that, though longer links strengthen a single-layer network, they make a multi-layer network more vulnerable. We further find that when ζ is longer than a certain critical value, , abrupt, discontinuous transitions take place, while for  the transition is continuous, indicating that the risk of abrupt collapse can be eliminated if the typical link length is shorter than .}
}
```