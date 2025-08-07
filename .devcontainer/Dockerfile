FROM mcr.microsoft.com/devcontainers/cpp:1-ubuntu-24.04

# Install cmake
RUN mkdir -p /usr/local
RUN cd /usr/local &&\
    wget https://github.com/Kitware/CMake/releases/download/v4.0.2/cmake-4.0.2-linux-x86_64.sh &&\
    chmod +x cmake-4.0.2-linux-x86_64.sh &&\
    ./cmake-4.0.2-linux-x86_64.sh --skip-license
    
# [Optional] Uncomment this section to install additional packages.
# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>
