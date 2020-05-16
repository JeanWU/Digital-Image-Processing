#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "string.h"
#define _USE_MATH_DEFINES
#include "math.h"

 
#define Size 256
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
        unsigned char gradient_map[Size][Size];
        int gradient_map2[Size][Size]; //for 2nd order 
        unsigned char edge_map[Size][Size];
	public:
		Image();

        double calcuate3point(int j, int k){
            int GR = 0, GC = 0;
            if(k > 0){
                GR = ImageData[j][k] - ImageData[j][k-1];
            }
            if(j < Size - 1){
                GC = ImageData[j][k] - ImageData[j+1][k];
            }
            return sqrt(pow(GR,2) + pow(GC,2));
        }

        double calcuate4point(int j, int k){
            int GR = 0, GC = 0;
            if((k < Size - 1) && (j < Size - 1)){
                GR = ImageData[j][k] - ImageData[j+1][k+1];
                GC = ImageData[j][k+1] - ImageData[j+1][k];
            }
            return sqrt(pow(GR,2) + pow(GC,2));
        }

        double calcuate9point(int j, int k, int maskType){
            //maskType = 1 => Prewitt; maskType = 2 => Sobel
            double GR = 0, GC = 0;
            if(j < Size-1 && j > 0 && k < Size-1 && k > 0){
				GR = round((ImageData[j-1][k+1]+maskType*ImageData[j][k+1]+ImageData[j+1][k+1] - \
							ImageData[j-1][k-1]+maskType*ImageData[j][k-1]+ImageData[j+1][k+1]) / (maskType+2));
				GC = round((ImageData[j-1][k-1]+maskType*ImageData[j-1][k]+ImageData[j-1][k+1] - \
							ImageData[j+1][k-1]+maskType*ImageData[j+1][k]+ImageData[j+1][k+1]) / (maskType+2));
			}
            return sqrt(pow(GR,2) + pow(GC,2));
        }
        
        int calcuate4or8neighbor(int i, int j, int neighborType){
            //neighborType=4, use mask4; neighborType=8, use mask8
            double mask4[3][3] = {{0,-1,0}, {-1,4,-1}, {0,-1,0}};
            double mask8[3][3] = {{-1,-1,-1},{-1,8,-1},{-1,-1,-1}};
            int filter_size = 3;
            int length = (filter_size-1)/2;
            
            double convolution_val = 0;
            for(int k=0; k<filter_size; k++){
                if( (i-length+k) < 0 || (i-length+k) >= Size) continue;
                for(int l=0; l<filter_size; l++){
                    if( (j-length+l) < 0 || (j-length+l) >= Size) continue;
                    switch(neighborType){
                        case 4:
                            convolution_val += mask4[k][l] * ImageData[i-length+k][j-length+l];
                            break;
                        case 8:
                            convolution_val += mask8[k][l] * ImageData[i-length+k][j-length+l];
                            break;
                        default: 
                            cout << "wrong neighbor type!" << endl;
                    }
                }
            }
            return round(convolution_val/neighborType);
        }

        void computeGradientMap(int type){
            //type = 3, use calculate 3 point; type = 4, use calculate 3 point;
            //type = 1, use Prewitt; type = 2, use Sobel
            //type =7, use calcuate4or8neighbor 4 neighbor
            //type =8, use calcuate4or8neighbor 8 neighbor
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    switch(type){
                        case 1: 
                            gradient_map[i][j] = calcuate9point(i,j,1);
                            break;
                        case 2:
                            gradient_map[i][j] = calcuate9point(i,j,2);
                            break;
                        case 3:
                            gradient_map[i][j] = calcuate3point(i,j);
                            break;
                        case 4:
                            gradient_map[i][j] = calcuate4point(i,j);
                            break;
                        case 7:
                            gradient_map2[i][j] = calcuate4or8neighbor(i,j,4);
                            break;
                        case 8:
                            gradient_map2[i][j] = calcuate4or8neighbor(i,j,8);
                            break;
                        default:
                            cout << "wrong type for calculating gradient" << endl;
                    }
                }
            }          
        }
        
        void out_hist(char out[100], int type){
            //count hist
            //type = 1, normal; type=2, for 2nd order
            int cnt[Size] = {0};
            int cnt2ndOrder[Size] ={0};
			for(int i = 0; i < Size; i++){
                for(int j = 0; j < Size; j++){
                    switch(type){
                        case 1:
                            cnt[gradient_map[i][j]]++;
                            break;
                        case 2:
                            if(gradient_map2[i][j]>=-128 && gradient_map2[i][j] <=127){
                                cnt2ndOrder[int(gradient_map2[i][j])+128]++;
                            }
                            break;
                        default:
                            cout << "wrong histogram type -1 !" << endl;
                    }
                }
            }
            //write csv file           
            ofstream myfile (out);
            if (myfile.is_open())
            {
                myfile << "value,count" << endl;
                switch(type){
                    case 1:
                        for(int i=0; i<Size; i++){
                            myfile << i << "," << cnt[i] << endl;
                        }
                        break;
                    case 2:
                        for(int i=-128; i<128; i++){
                            myfile << i << "," << cnt2ndOrder[i+128] << endl;
                        }
                        break;
                    default:
                        cout << "wrong histogram type -2 !" << endl;
                }
                myfile.close();
            }
            else {
                cout << "Unable to open file";
                exit(1);
            }
        }
        
        void generateEdgeMap(int threshold, int type){
            switch(type){
                case 1:
                    for(int i=0; i<Size; i++){
                        for(int j=0; j<Size; j++){
                            if(gradient_map[i][j] >= threshold){
                                edge_map[i][j] = 255;
                            }
                        }
                    }
                    break;
                case 2:
                    for(int i=0; i<Size; i++){
                        for(int j=0; j<Size; j++){
                            if(abs(gradient_map2[i][j]) <= threshold) gradient_map2[i][j] = 0;
                        }
                    }
                    for(int i=0; i<Size; i++){
                        for(int j=0; j<Size; j++){
                            if(gradient_map2[i][j] == 0){
                                if(gradient_map2[i-1][j]*gradient_map2[i+1][j] < 0) edge_map[i][j] = 255;
                                if(gradient_map2[i-1][j-1]*gradient_map2[i+1][j+1] < 0) edge_map[i][j] = 255;
                                if(gradient_map2[i-1][j+1]*gradient_map2[i+1][j-1] < 0) edge_map[i][j] = 255;
                                if(gradient_map2[i][j-1]*gradient_map2[i][j+1] < 0) edge_map[i][j] = 255;
                            }
                        }
                    }
                    break;
            }
            memcpy(ImageData, edge_map, sizeof(unsigned char)*Size*Size);
        }

        int Connect(int x, int y){
            if(x < 0 || x > Size-1 || y < 0 || y > Size-1)  return 0;
            for(int i = -1; i < 2; i++)
                for(int j = -1; j < 2; j++){
                    if(x+i < 0 || x+i > Size-1 || y+j < 0 || y+j > Size-1)  continue;
                    if(edge_map[x+i][y+j] == 2){
                        edge_map[x][y] = 2;
                        ImageData[x][y] = 255;
                        return 2;
                    }
                }
            for(int i = -1; i < 2; i++)
                for(int j = -1; j < 2; j++){
                    if(x+i < 0 || x+i > Size-1 || y+j < 0 || y+j > Size-1)  continue;
                    if(edge_map[x+i][y+j] == 1){
                        edge_map[x+i][y+j] = Connect(x+i, y+j);
                    }
                }
        }

        void Canny(double TL, double TH){
            int filter_sz = 5;
			int len = (filter_sz-1)/2;
			unsigned char tmp[Size][Size];
			
			double gaussianFilter[5][5] = {{2,4,5,4,2},{4,9,12,9,4},{5,12,15,12,5},{4,9,12,9,4},{2,4,5,4,2}};
			for(int i = 0; i < Size; i++){
				for(int j = 0; j < Size; j++){
					double sum = 0, cnt = 0;
					for(int k = 0; k < filter_sz; k++){
						if(i-len+k < 0 || i-len+k > Size-1)	continue;
						for(int l = 0; l < filter_sz; l++){
							if(j-len+l < 0 || j-len+l > Size-1)	continue;
							sum += int(ImageData[i][j])*gaussianFilter[k][l];
							cnt += gaussianFilter[k][l];
						}
					}
					tmp[i][j] = round(sum/cnt);
				}
			}
			memcpy(ImageData, tmp, sizeof(unsigned char)*Size*Size);
			int grad_map[Size][Size] = {{0}};
			float theta[Size][Size] = {{0}};
			for(int i = 0; i < Size; i++){
                for(int j = 0; j < Size; j++){
					if(i == Size-1 || j == 0)	continue;
					double GR = 0, GC = 0;
					GC = ImageData[i][j] - ImageData[i+1][j];
					GR = ImageData[i][j] - ImageData[i][j-1];
					grad_map[i][j] = sqrt(pow(GC,2)+pow(GR,2));
                    theta[i][j] = atan2(GC, GR) * 180 / M_PI;
				}
            }
				
			int edge_map[Size][Size] = {{0}};
            for(int i = 1; i < Size-1; i++){
                for(int j = 1; j < Size-1; j++){
                    if(theta[i][j] < 22.5 && theta[i][j] > -22.5){
                        if(grad_map[i][j] > grad_map[i][j-1] && grad_map[i][j] > grad_map[i][j+1])
                            edge_map[i][j] = grad_map[i][j];
                    }
                    else if(theta[i][j] > 22.5 && theta[i][j] < 67.5){
                        if(grad_map[i][j] > grad_map[i-1][j+1] && grad_map[i][j] > grad_map[i+1][j-1])
                            edge_map[i][j] = grad_map[i][j];
                    }
                    else if(theta[i][j] < -22.5 && theta[i][j] > -67.5){
                        if(grad_map[i][j] > grad_map[i+1][j+1] && grad_map[i][j] > grad_map[i-1][j-1])
                            edge_map[i][j] = grad_map[i][j];
                    }
                    else{
                        if(grad_map[i][j] > grad_map[i-1][j] && grad_map[i][j] > grad_map[i+1][j])
                            edge_map[i][j] = grad_map[i][j];
                    }
                }
            }
            for(int i = 0; i < Size; i++)
                for(int j = 0; j < Size; j++){
                    if(edge_map[i][j] >= TH)                            edge_map[i][j] = 2;
                    else if(edge_map[i][j] < TH && edge_map[i][j] > TL) edge_map[i][j] = 1;
                    else                                                edge_map[i][j] = 0;
                }
            memset(ImageData, 0, Size*Size*sizeof(unsigned char));
            for(int i = 1; i < Size-1; i++)
                for(int j = 1; j < Size-1; j++){
                    if(edge_map[i][j] == 1) Connect(i, j);
                    if(edge_map[i][j] == 2) ImageData[i][j] = 255;
                }
		}
		
		void output(char out[100]){
			FILE *fout;
			if (!(fout=fopen(out,"wb"))){
				cout<<"Cannot open file!"<<endl;
				exit(1);
			}
			fwrite(ImageData, sizeof(unsigned char), Size*Size, fout);
			fclose(fout);
		}
        
};
Image::Image(){
	memcpy(ImageData, Original_img, sizeof(unsigned char)*Size*Size);
	memset(gradient_map, 0, sizeof(unsigned char)*Size*Size);
    memset(edge_map, 0, sizeof(unsigned char)*Size*Size);
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
    
    Image sample1E1, sample1E2, sample1E3;

    sample1E1.computeGradientMap(3);
    sample1E1.generateEdgeMap(40,1);
    sample1E1.output(argv[2]);
    sample1E2.computeGradientMap(8);
    sample1E2.generateEdgeMap(15,2);
    sample1E2.output(argv[3]);
    sample1E3.Canny(10, 40);
    sample1E3.output(argv[4]);


    if (!(fin=fopen(argv[5],"rb"))){
		cout<<"Cannot open file!"<<endl;
		exit(1);
	}
	fread(Original_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);
    
    Image sample2E1, sample2E2, sample2E3;

    sample2E1.computeGradientMap(3);
    sample2E1.generateEdgeMap(80,1);
    sample2E1.output(argv[6]);
    sample2E2.computeGradientMap(8);
    sample2E2.generateEdgeMap(30,2);
    sample2E2.output(argv[7]);
    sample2E3.Canny(50, 70);
    sample2E3.output(argv[8]);

    if (!(fin=fopen(argv[9],"rb"))){
		cout<<"Cannot open file!"<<endl;
		exit(1);
	}
	fread(Original_img, sizeof(unsigned char), Size*Size, fin);
	fclose(fin);
    
    Image sample3E1, sample3E2, sample3E3;

    sample3E1.computeGradientMap(3);
    sample3E1.generateEdgeMap(40,1);
    sample3E1.output(argv[10]);
    sample3E2.computeGradientMap(8);
    sample3E2.generateEdgeMap(15,2);
    sample3E2.output(argv[11]);
    sample3E3.Canny(30, 60);
    sample3E3.output(argv[12]);
    
	return 0;
}
