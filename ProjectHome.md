An image processing algorithm for labeling regions in an open-cell cellular structure. We have currently published the source code to the primary tool (DistLabel.java) and accompanying files. It can currently be used as part of another Java program which can produce the distance map and object masks as a 3D array represented in a linear fashion
```
 [x][y][z] -> [(z*width+y)*height+x] 
```
.
## Upcoming ##
We have migrated the project including latest version of the code to https://bitbucket.org/skicavs/tipl-public This includes much more code to read, view, and save images as well as the most up to date version of the tools.