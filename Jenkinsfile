pipeline {
    //agent any
    agent {
        docker {
                image 'ubuntu:18.04' // 'qnib/pytest'
        }
    }

    environment {
        CODE_COV_TOKEN = credentials('CODE_COV_TOKEN')
    }

    stages {
        stage('Build') {
            steps {
                sh script: '''apt update && apt install  -y python3.8 apt-utils python3-pip lsb-release wget curl software-properties-common && \
                 bash -c "$(wget -O - https://apt.llvm.org/llvm.sh)" && pip3 install --upgrade pip && pip3 install pytest && pip3 install pytest-cov && \
                 pip3 install codecov && pip3 install numpy  && pip3 install pytest-cov && pip3 install -r requirements.txt && python3 setup.py install'''
            }
        }
        stage('Test') {
            steps {
                //withCredentials([file(credentialsId: 'secret', variable: 'CODE_COV_TOKEN')]) {
                withCredentials([string(credentialsId: 'CODE_COV_TOKEN', variable: 'CODE_COV_TOKEN')]) {
                    sh script: '''
                    python3 -m pytest tests/ --cov=./norec4dna && \
                      curl -s https://codecov.io/bash | bash -s - -t $CODE_COV_TOKEN'''
                }
            }
        }
        stage('Deploy') {
            steps {
                echo 'Deploying....'
            }
        }
    }
}