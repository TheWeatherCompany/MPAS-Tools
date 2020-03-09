pipeline {
    environment {
        SOURCE_DIR='/gpfs/data1/storage/home/regression/.jenkins/workspace/graf_repos_src'
    }
    agent any
    stages {
        stage('Make') {
            steps {
                sh 'pwd'
                sh 'cd ${SOURCE_DIR}/GlobalNWP; pwd; source ./dtrc dyeus; cd ../MPAS-Tools; cd limited_area/mpas_to_mpas/; ./build_wsc.sh'
            }
        }
    }
}
