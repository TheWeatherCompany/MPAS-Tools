pipeline {
    agent any
    stages {
        stage('Make') {
            steps {
                sh 'pwd'
                sh 'cd ../graf_repos_src/GlobalNWP; pwd; source ./dtrc dyeus; cd ../MPAS-Tools; cd limited_area/mpas_to_mpas/; ./build_wsc.sh'
            }
        }
    }
}
