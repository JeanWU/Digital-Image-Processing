#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "string.h"
#include <cmath>
 
#define Size 512
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
        unsigned char ImageData_bk[Size][Size];
        double *** M;
        double *** T;
        double ** kmeansCentroids;
        unsigned int ** kmeansCluster;
        unsigned int ** kmeansClusterBefore;
        int CLUSTERCOUNT = 4;
        unsigned char ** texture0;
        unsigned char ** texture1;
        unsigned char ** texture2;
        unsigned char ** texture3;
        //unsigned char ** outputImage;

	public:
		Image();

        void laws_method(int windowSize){
            M = allocMemory3DDouble(M,Size,Size,9);
            T = allocMemory3DDouble(M,Size,Size,9);
            
            for(int x=1; x<Size-1; x++){
                for(int y=1; y<Size-1; y++){
                    for(int k=0; k<9; k++){
                        M[x][y][k] = convolution(x,y,k);
                    }
                }
            }

            signed int windowStart = - int ((windowSize -1)/2);
            signed int windowEnd = int ((windowSize -1)/2);
            
            for(int x=0; x<Size-1; x++){
                for(int y=0; y<Size-1; y++){
                    for(int k=0; k<9; k++){
                        for(signed int m=windowStart; m<=windowEnd; m++){
                            for(signed int n=windowStart; n<=windowEnd; n++){
                                signed int rowIndex = x + m;
                                signed int colIndex = y + n;
                                if(rowIndex < 0) {
                                    rowIndex = 0;
                                }
                                if(colIndex < 0) {
                                    colIndex = 0;
                                }
                                if(rowIndex > Size-1) {
                                    rowIndex = Size-1;
                                }
                                if(colIndex > Size-1) {
                                    colIndex = Size-1;
                                }
                                T[x][y][k] += pow(abs(M[rowIndex][colIndex][k]),2);
                            }
                        }
                        T[x][y][k] = T[x][y][k] / (windowSize*windowSize);
                        //T[x][y][k] = sqrt(T[x][y][k]);
                    }
                }
            }
        }

        float convolution(int x, int y, int filterType){
            float H[9][3][3] = {{{1/36,2/36,1/36},{2/36,4/36,2/36},{1/36,2/36,1/36}},
                                {{1/12,0,-1/12},{2/12,0,-2/12},{1/12,0,-1/12}},
                                {{-1/12,2/12,-1/12},{-2/12,4/12,-2/12},{-1/12,2/12,-1/12}},
                                {{-1/12,-2/12,-1/12},{0,0,0},{1/12,2/12,1/12}},
                                {{1/4,0,-1/4},{0,0,0},{-1/4,0,1/4}},
                                {{-1/4,2/4,-1/4},{0,0,0},{1/4,-2/4,1/4}},
                                {{-1/12,-2/12,-1/12},{2/12,4/12,2/12},{-1/12,-2/12,-1/12}},
                                {{-1/4,0,1/4},{2/4,0,-2/4},{-1/4,0,1/4}},
                                {{1/4,-2/4,1/4},{-2/4,4/4,-2/4},{1/4,-2/4,1/4}}};
            float sum = 0;
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    sum += ImageData[x+i-1][y+j-1]*H[filterType][i][j];
                }
            }
            return sum;
        }
        
        void initializeKmeansCentroids(){
            kmeansCentroids = allocMemory2DDouble(kmeansCentroids, CLUSTERCOUNT, 9);
            for(int i=0; i< CLUSTERCOUNT; i++) {
                for(int j=0; j<9; j++) {
                    if(i==0) {
                        kmeansCentroids[i][j] = T[110][100][j];
                    }
                    else if(i==1) {
                        kmeansCentroids[i][j] = T[80][400][j];
                    }
                    else if(i==2) {
                        kmeansCentroids[i][j] = T[320][210][j];

                    }
                    else if(i==3) {
                        kmeansCentroids[i][j] = T[400][430][j];
                    }
                }
            }
        }

        void kmeansClustering(){
            kmeansCluster = allocMemory2DInt(kmeansCluster, Size, Size);
            kmeansClusterBefore = allocMemory2DInt(kmeansClusterBefore, Size, Size);
            bool flag = true;
            int iteration = 0;
            //for(int iterationNO=0; iterationNO < ITERATIONCOUNT; iterationNO++) {
            while(flag){
                //cout << "iteration :" << iterationNO << endl;
                for(int x=0; x<Size; x++){
                    for(int y=0; y<Size; y++){
                        double * tempDistance = new double [CLUSTERCOUNT]();

                        //calculate distance from each point to all the cluster centroids
                        for (int clusterNo=0; clusterNo < CLUSTERCOUNT; clusterNo++) {
                            long double tempSum = 0;
                            for (int i = 0; i < 9; i ++){
                                tempSum = tempSum + pow((T[x][y][i] - kmeansCentroids[clusterNo][i]) , 2);
                            }
                            tempSum = sqrt(tempSum);
                            tempDistance[clusterNo] = tempSum;
                        }

                        //get minimun distance and its corresponding centroids
                        long double minDist = tempDistance[0];
                        unsigned int minDistCluster = 0;
                        for(int i=1; i< CLUSTERCOUNT; i++) {
                            if( tempDistance[i] < minDist ) {
                                minDist = tempDistance[i];
                                minDistCluster = i;
                            }
                        }
                        kmeansCluster[x][y] = minDistCluster;
                    }
                }              
                
                //Update Centroids
                for(int clusterNo=0; clusterNo < CLUSTERCOUNT; clusterNo++) {
                    
                    int updateCount = 0;
                    for(int x=0; x<Size; x++){
                        for(int y=0; y<Size; y++){
                            if( kmeansCluster[x][y] == clusterNo) {
                                if( updateCount == 0) {
                                    for(int i=0; i < 9; i++) {
                                        kmeansCentroids[clusterNo][i] = T[x][y][i];
                                        //cout << kmeansCentroids[clusterNo][i] << endl;
                                    }
                                }
                                else {
                                    for(int i=0; i < 9; i++) {
                                        kmeansCentroids[clusterNo][i] = kmeansCentroids[clusterNo][i] + T[x][y][i];
                                        //cout << kmeansCentroids[clusterNo][i] << endl;
                                    }
                                }
                                updateCount++;
                            }
                        }
                    }
                    
                    for(int i=0; i < 9; i++) {
                        if(updateCount > 0){
                            kmeansCentroids[clusterNo][i] = kmeansCentroids[clusterNo][i] / ((double) (updateCount));
                        }
                        //cout << kmeansCentroids[clusterNo][i] << endl;
                    }
                }

                //compare kmeansCluster and kmeansClusterBefore; the iteration converges if equals
                bool isDifference = false;
                for(int i=0; i<Size; i++){
                    for(int j=0; j<Size; j++){
                        if(kmeansClusterBefore[i][j] != kmeansCluster[i][j]){
                            isDifference = true;
                            continue;
                        }
                    }
                }
                memcpy(kmeansClusterBefore, kmeansCluster, sizeof(double)*Size*Size);
                if(!isDifference)   flag = false;

                iteration = iteration + 1;
                // if(iteration == 20){
                //     flag = false;
                // }
            }
            //cout << "iteration: " << iteration << endl;
            freeMemory3DDouble(T,Size,Size,9);
            freeMemory2DDouble(kmeansCentroids,4,9);
            freeMemory2DInt(kmeansClusterBefore,Size,Size);
        }

        void clustersToGraylevels() {
            unsigned char graylevel[4] = {0,64,128,255};
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    //cout << kmeansCluster[tempCount] << endl;
                    for(int k=0; k < CLUSTERCOUNT; k++) {
                        if(kmeansCluster[i][j] == k) {
                            ImageData[i][j] = graylevel[k];
                            //cout << kmeansCluster[tempCount] << ' ' << (int)graylevel[k] << endl;
                        }
                    }
                }
            }
        }

        void swap_feature(){
            texture0 = allocMemory2D(texture0,50,50);
            texture1 = allocMemory2D(texture1,50,50);
            texture2 = allocMemory2D(texture2,50,50);
            texture3 = allocMemory2D(texture3,50,50);
            //outputImage = allocMemory2D(outputImage,Size,Size);
            unsigned char outputImage[Size][Size] = {{0}};

            
            for(int x=0; x<50; x++){
                for(int y=0; y<50; y++){
                    texture0[x][y] = ImageData_bk[x][y];
                    texture1[x][y] = ImageData_bk[x][400+y];
                    texture2[x][y] = ImageData_bk[420+x][y];
                    texture3[x][y] = ImageData_bk[400+x][400+y];
                }
            }

            for(int x=0; x<Size; x+=50){
                for(int y=0; y<Size; y+=50){
                    for(int i=0; i<50; i++){
                        for(int j=0; j<50; j++){
                            if((x+i>Size-1) || (y+j>Size-1))    continue;
                            if(kmeansCluster[x+i][y+j] == 0)    outputImage[x+i][y+j] = texture1[i][j];
                            if(kmeansCluster[x+i][y+j] == 1)    outputImage[x+i][y+j] = texture2[i][j];
                            if(kmeansCluster[x+i][y+j] == 2)    outputImage[x+i][y+j] = texture3[i][j];
                            if(kmeansCluster[x+i][y+j] == 3)    outputImage[x+i][y+j] = texture0[i][j];                       
                        }
                    }
                }
            }
            memcpy(ImageData, outputImage, sizeof(unsigned char)*Size*Size);
        }
        		
		void output(char out[100]){
			FILE *fout;
			if (!(fout=fopen(out,"wb"))){
				cout<<"Cannot open file!"<<endl;
				exit(1);
			}
			fwrite(ImageData, sizeof(unsigned char), Size*Size, fout);
            //fwrite(outputImage, sizeof(unsigned char), Size*Size, fout);
			fclose(fout);
		}

        void output_feature_vector(){
            for(int i=0; i<9; i++){
                char outputFileName[10];
                sprintf(outputFileName, "M_%d.raw", i+1);
                FILE *fout;
                if (!(fout=fopen(outputFileName,"wb"))){
                    cout<<"Cannot open file!"<<endl;
				    exit(1);
                }
                fwrite(M[i], sizeof(unsigned char), Size*Size*1, fout);
                fclose(fout);
            }

            for(int i=0; i<9; i++){
                char outputFileName[10];
                sprintf(outputFileName, "T_%d.raw", i+1);
                FILE *fout;
                if (!(fout=fopen(outputFileName,"wb"))){
                    cout<<"Cannot open file!"<<endl;
				    exit(1);
                }
                fwrite(T[i], sizeof(unsigned char), Size*Size*1, fout);
                fclose(fout);
            }
            freeMemory3DDouble(M,Size,Size,9);
        }

        double *** allocMemory3DDouble(double *** image3D, int ROW, int COL, int BYTESPERPIXEL) {
            image3D = new double **[ROW]();
            for(int i=0; i < ROW; i++) {
                image3D[i] = new double *[COL]();
                for(int j=0; j < COL; j++) {
                    image3D[i][j] = new double [BYTESPERPIXEL]();
                    for(int k=0; k < BYTESPERPIXEL; k++) {
                        image3D[i][j][k] = 0;
                    }
                }
            }
            return image3D;
        }

        double ** allocMemory2DDouble(double ** image2D, int ROW, int COL) {
            image2D = new double *[ROW]();
            for(int i=0; i < ROW; i++) {
                image2D[i] = new double [COL]();
                for(int j=0; j < COL; j++) {
                    image2D[i][j] = 0;
                }
            }
            return image2D;
        }
        
        unsigned int ** allocMemory2DInt(unsigned int ** image2D, int row, int col){
            image2D = new unsigned int *[row]();
            long int temp_count = 0;
            for(int i=0; i < row; i++) {
                image2D[i] = new unsigned int [col]();
                for(int j=0; j < col; j++) {
                    image2D[i][j] = 0;
                    temp_count = temp_count + 1;
                }
            }
            return image2D;
        }

        unsigned char ** allocMemory2D(unsigned char ** image2D, int row, int col){
            image2D = new unsigned char *[row]();
            long int temp_count = 0;
            for(int i=0; i < row; i++) {
                image2D[i] = new unsigned char [col]();
                for(int j=0; j < col; j++) {
                    image2D[i][j] = 0;
                    temp_count = temp_count + 1;
                }
            }
            return image2D;
        }

        void freeMemory3DDouble(double *** image3D, int row, int col, int bytesPerPixel) {
            for (int i=0; i<row; i++) {
                for (int j=0; j<col; j++) {
                    delete[] image3D[i][j];
                }
                delete[] image3D[i];
            }
            delete[] image3D;
        }

        void freeMemory2DDouble(double ** image2D, int row, int col) {
            for (int i=0; i<row; i++) {
                delete[] image2D[i];
            }
            delete[] image2D;
        }

        void freeMemory2DInt(unsigned int ** image2D, int row, int col) {
            for (int i=0; i<row; i++) {
                delete[] image2D[i];
            }
            delete[] image2D;
        }
};
Image::Image(){
	memcpy(ImageData, Original_img, sizeof(unsigned char)*Size*Size);
    memcpy(ImageData_bk, Original_img, sizeof(unsigned char)*Size*Size);
}

int main(int argc, char **argv)
{
	FILE *fin;
	if (!(fin=fopen(argv[1],"rb"))){
		cout<<"Cannot open file!"<<endl;
		exit(1);
	}
	fread(Original_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);
    
    Image featureImage;
    featureImage.laws_method(51);
    featureImage.output_feature_vector();
    featureImage.initializeKmeansCentroids();
    featureImage.kmeansClustering();
    featureImage.clustersToGraylevels();
    featureImage.output(argv[2]);
    featureImage.swap_feature();
    featureImage.output(argv[3]);

	return 0;
}