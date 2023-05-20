# SVanalysis_STARProtocols
A protocol for applying low-coverage whole-genome sequencing data in structural variation studies


For complete details on the use and execution of this protocol, please refer to Liu et al. (2023)[1]


Reference


[1] Liu, Q., Yang, K., Xie, B., Gao, Y., Xu, S., and Lu, Y. (2023). Mapping Structural Variations in Haemaphysalis longicornis and Rhipicephalus microplus Reveals Vector–Pathogen Adaptation. iScience 26, 106398

# How to install docker image

1. build image: docker build -t svanalysis_starprotocols:1.0.0 ./ --no-cache
2. make a new container from above image: docker run -it --name test -p 6770:22 svanalysis_starprotocols:1.0.0 /bin/bash
3. login through ssh: ssh root@0.0.0.0 -p 6770
4. you can do anything like in Linux server
