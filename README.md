# Augmented-Reality

In this project, a simple augmented reality application has been developed. The deliverable is a video that contains several virtual object models as if they exist in the real
world. Furthermore, you will be able to specify pixel positions to place an arbitrary
object.
A video with AprilTags (https://april.eecs.umich.edu/software/apriltag) in each frame is our input. These tags are usually used in robotics for determining the pose of the
camera. 4 corners coordinates (in pixel) as well as the size of the tags are also provided. 

Camera poses have then been recovered with two different approaches: 
1) solving the Perspective-N-Point (PnP) problem with coplanar assumption
2) solving the Persepective-three-point (P3P) and the Procrustes problem.

After retrieving the 3D relationship between the camera and world, we can place an
arbitrary objects in the scene.

![image](https://drive.google.com/uc?export=view&id=12RKcF5_6rjIMhtI6V4MtwTPQtAdVa7LG)
![image](https://drive.google.com/uc?export=view&id=GdhrDuaCke1QHsdIZL6H1airMs1pjf)
