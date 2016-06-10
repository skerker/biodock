BIODOCK
====
A Docker container that runs bioinformatics software.
Biodock was created for bioinformatic analysis at the Hinton Lab.


General updates
====

    11.05.2016:   Git created.
    03.06.2016:   biodock updated - GitLab and DockerHub projects are now linked, the Docker image now runs using ubuntu and contains a completely new list of software

-----------




Contents
====

1.	Running biodock
---

  If you are on a Mac, start your Docker machine:

  `docker-machine start MACHINE-NAME`

  `eval "$(docker-machine env MACHINE-NAME)"`


  Pull the biodock image from Docker Hub:

  `docker pull wpmr/biodock:latest`


  Alternatively, clone this git and build the biodock image from the Dockerfile:

  `git clone https://gitlab.com/will_rowe/biodock.git`

  `cd biodock`

  `docker build -t wpmr/biodock:latest .`


  Launch the Docker container, making sure to mount a volume (allowing you to transfer data in and out of the container):

  `docker run -itP -m 8g --name biodock -v /some/path/on/host/:/MOUNTED_VOLUME wpmr/biodock:latest`

  + -i = keep STDIN open even if not attached

  + -t = allocate a pseudo-tty

  + -P = publish all exposed ports to the host interfaces

  + -m = memory limit (8gb)

  + --name = name for container at runtime (easy to use for later exec commands)

  + -v = bind mount a volume (for data transfer etc. between container and host machine). Usage-> [host-src:]container-dest[:<options>]. The comma-delimited `options` are [rw|ro], [z|Z], [[r]shared|[r]slave|[r]private], and [nocopy].




2.	Usage
----

  This Docker container (and git repo) is intended to provide a standardised environment for running bioinformatics pipelines, scripts and software.


  The container will launch bash by default, all software is in the path and scripts from the git repo are in /opt/SCRIPT_bin (also in path)


  A few helpful commands for managing the container:

  + Once exited, you can re-enter the container using the exec command:

  `docker exec -it [CONTAINER ID] bash`

  + View all containers (both running and stopped) using:

  `docker ps -a`

  + Stop or remove all containers

  `docker stop $(docker ps -a)`
  
  `docker rm $(docker ps -a)`




3. Notes
----

  + The CPUs available to Docker are limited by the host machine running docker, so set the virtual machine to have the required number before running Docker.

  + The Kernel scheduler will handle the resource contention in the case of multiple containers requiring multiple cores.
