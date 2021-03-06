# dynamics-simulations-of-human-movement
OVERVIEW OF THR PROJECT

the project aims to understand the fundamental dynamics of human movement, and proposes a model to capture several important features of human movement.

The mechanisms of human movement are very complicated, there is no general rule of thumb that explains what's happening in real work. However, based on the real-observed data in preschool classrooms, we do find some general patterns that can be helpful to understand the mechanism and build an explainable model.

      1.distances between children in the classroom are not randomly distributed. Data shows that children tend to interact with others at a "comfortable" distance, and radial distribution function g(r) indicates a peak around this particular distance.

      2.relative orientations of children and their partners are not randomly distributed. Data also shows that children tend to be oriented in a "shoulder-to-shoulder" way when the distance between is short. Whereas children are more likely to be oriented in a "face-to-face" way at a relatively larger distance.

Based on the above-observed patterns, we build a model that assumes children behave as ferromagnets at a short distance, while antiferromagnets at a relatively larger distance. The results show this model well captures several feathers of human movement.


REAL DATA VISUALIZATION

Radio Frequency Identification technology (RFID)vallows for efficient capture of each child’s location and movement throughout the classroom and informs our understanding of which children tend to be in proximity of each other.        
      
<img src="findings/traj-300sec.png" width=600 >\
<em>Figure 1. continuous measurement of real human locations and orientations </em>

Figure 1 shows the trajectories of the nine children for 300 seconds during one period of class time. Using tags worn by children, we can also find orientation for each kid data collection. The arrow represents the orientation of the child. We find that the children are not uniformly walking/running the whole classroom space, they wander in a small area for a certain amount of time when they are interacting with other children in social contact, the motions are trapped by localization.  For instance, once children were approaching others (see several clusters in this center of Fig. 1, each cluster shows a trajectory of a child),  children kept roaming around for a while.


