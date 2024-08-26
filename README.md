# Data-driven Force Observer for Human-Robot Interaction using Gaussian Processes
This repository contains an implementation of the paper [1]:

[1] S. Tesfazgi*, M. Keßler*, E. Trigili, A. Lederer and S. Hirche. "Data-driven Force Observer for Human-Robot Interaction with Series Elastic Actuators using Gaussian Processes" in IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2024 [arxiv](https://arxiv.org/abs/2405.08711)

## Repository structue
The "simulation files" directory contains implementations of the GP-AKF method for estimation of external torques in simulations. The simulation is implemented for a 1-DoF and a 2-DoF SEA exoskeleton. Human passive dynamics can be simulated as simple additive inertia, mass and damping models. Finally, an implementation of the GP-AKF with the guaranteed error bounds from the paper is given for the 2-DoF simulation.

## Code Structure
The source code of the simulation uses two generic classes "GP-AEKF.m" and "robot_system.m". These are used identically in each simulation implementation. The classes are fully documented and information on usage, parameters, inputs and outputs of functions is given in the respective file. 

The first class called "GP-AEKF" implements the GP-augmented torque estimation method from the paper. For the implementation, mixture of experts LoG-GPs and the extended Kalman filter are used. The most important functions are:  
-"train_loggp": Update LoG-GP model with a new sample. Load side angle/velocity/acceleration ($q,\dot{q},\ddot{q}$) and motor torques ($\tau_m$) must be given.  
-"initialize": Initializes Kalman filter with initial states $x_0$ and covariance $P_0$. Must be called before estimation.  
-"predict": Performs Kalman filter prediction based on enhanced-GP model. The motor torque input $tau_m$ must be given. The predicted state and covariance are saved as private properties of the object.  
-"update": Performs filtering step of the KF based on measurements. The measurements $y$ must be given. The estimate and estimate covariance can be retrieved from the properties $\hat{x}$ and $P$, where they are saved after calling "update".  

The second class called "robot_system.m" implements the system dynamics based on given parameters. This provides functions to compute inverse dynamics and forward dynamics of the nominal system based on symbolic expressions. The symbolic expressions are converted to much faster Matlab functions based on the given parameters. Examples of how the parameters (properties of the class) can be set are given in the respective "config.m" files of each application directory. The "robot_system.m" class is required and used by the "GP-AEKF.m" class. Therefore, an instantiated object of "robot_system.m" must be given to the property "nom_model" of the "GP-AEKF.m" class. Important functions are:
- "fwd_dynamics" gives the nominal forward dynamics used by the GP-AEKF predictions. It accepts motor torque inputs tau_m
- "inverse dynamic" gives the nominal inverse dynamics of the load side. These are used to calculate the torque residual learned by the LoG-GP in the "GP-AEKF.m" class.

**Important:**
The class "robot_system.m" must be updated for different applications. At the moment, it only contains symbolic expressions for 1 and 2-DoF robots with a cylindrical intertia model. For different DoFs or models, the expressions must be changed or added! In more detail, the load side interia matrix **$M(q)$** must be inserted in "get_M_robot_sym" as a symbolic Matlab expression in terms of the generalized coordinates **q** and the classes parameters. Furthermore, the expressions for the Coriolis/Centrifugal forces matrix **$C(q,\dot{q})$** and the gravity torque **$G(q)$** must be set as symbolic Matlab expressions in terms of the generalized coordinates **$q,\dot{q}$** and the classes parameters.

---

If you found this software useful for your research, consider citing us.
```
@misc{tesfazgi2024,
      title={Data-driven Force Observer for Human-Robot Interaction with Series Elastic Actuators using Gaussian Processes}, 
      author={Samuel Tesfazgi and Markus Keßler and Emilio Trigili and Armin Lederer and Sandra Hirche},
      year={2024},
      eprint={2405.08711},
      archivePrefix={arXiv},
      primaryClass={cs.RO},
      url={https://arxiv.org/abs/2405.08711}, 
}
```
