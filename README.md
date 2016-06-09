BioDock git
====


General updates
====

    11.05.2016:   Git created.
    03.06.2016:   BioDock completely updated - GitLab and DockerHub projects are now linked, the Docker image now contains a completely new list of software

-----------


Contents
====

1.	Running BioDock
---

+ Clone this git and pull the BioDock image from Docker Hub:

`docker pull wpmr/biodock:latest`

+ Launch the Docker container, making sure to mount a volume:

`docker run -itP -m 8g --name BioDock -v /Users/willrowe/OneDrive/POSTDOC/0004_bioinformatics/dockyard/:/DIR_MOUNT wpmr/biodock:latest bash`

+ -i = keep STDIN open even if not attached

+ -t = allocate a pseudo-tty

+ -P = publish all exposed ports to the host interfaces

+ -m = memory limit (8gb)

+ --name = name for container at runtime (easy to use for later exec commands)

+ -v = bind mount a volume (for data transfer etc. between container and host machine). Usage-> [host-src:]container-dest[:<options>]. The comma-delimited `options` are [rw|ro], [z|Z], [[r]shared|[r]slave|[r]private], and [nocopy].


2.	Usage
----

+ This Docker container image and git is intended to house a standardised environment for running bioinformatics scripts and software.

+ The container will have your dockyard directory (containing 00_SCRIPTS, SCRATCH and this README) mounted at /DIR_MOUNT.

+ Once exited, you can re-enter the container using the exec command:

`docker exec -it [CONTAINER ID] bash`

+ View all containers (both running and stopped) using:

`docker ps -a`


3. Notes
----

+ The CPUs available to Docker are limited by the host machine running docker, so set the virtual machine to have the required number before running Docker.

+ The Kernel scheduler will handle the resource contention in the case of multiple containers requiring multiple cores.
