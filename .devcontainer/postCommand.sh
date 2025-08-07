set -ue

ENV_FILE=".devcontainer/environment.yml"
ENV_NAME=$(grep '^name:' ${ENV_FILE} | cut -d ' ' -f2)
echo Activating/creating env $ENV_NAME
conda env list 
conda env create -f ${ENV_FILE} || echo 'Environment already exists'
echo "conda activate $ENV_NAME" >> ~/.bashrc