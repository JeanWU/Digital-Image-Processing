#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "string.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <list>

#define Size 256
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
		unsigned char lowPassResult[Size][Size];
		unsigned char circle[Size][Size];
		unsigned char swirl[Size][Size];
		list<double> checkXPos;
		list<double> checkYPos;

	public:
		Image();

		void high_pass(int maskType){
			//maskType=1, use mask1; maskType=2, use mask2
			int filter_sz = 3;
			int len = (filter_sz-1)/2;
			int mask1[3][3] = {{0,-1,0},{-1,5,-1},{0,-1,0}};
			int mask2[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};
			int temp[Size][Size] = {{0}};
			for(int i = 0; i < Size; i++){
				for(int j = 0; j < Size; j++){
					double sum = 0;
					for(int k = 0; k < filter_sz; k++){
						if(i-len+k < 0 || i-len+k > Size-1)	continue;
						for(int l = 0; l < filter_sz; l++){
							if(j-len+l < 0 || j-len+l > Size-1)	continue;
							switch(maskType){
								case 1:
									sum += int(ImageData[i][j])*mask1[k][l];
									break;
								case 2:
									sum += int(ImageData[i][j])*mask2[k][l];
									break;
								default:
									cout << "wrong maskType" << endl;
							}							
						}
					}
					temp[i][j] = sum;
				}
			}
			for(int i=0; i<Size; i++){
				for(int j=0; j<Size; j++){
					ImageData[i][j] = temp[i][j];
				}
			}
		}
        
        void low_pass(double weight){
            int filter_size = 3;
            int length = (filter_size-1)/2;
            int temp[Size][Size];
            double mask[3][3] = {{1,weight,1}, {weight,pow(weight,2),weight}, {1,weight,1}};
                    
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    double convolution_val = 0;
                    for(int k=0; k<filter_size; k++){
                        if( (i-length+k) < 0 || (i-length+k) >= Size) continue;
                        for(int l=0; l<filter_size; l++){
                            if( (j-length+l) < 0 || (j-length+l) >= Size) continue;
                            convolution_val += mask[k][l] * ImageData[i-length+k][j-length+l];
                        }
                    }
                    temp[i][j] = convolution_val;
                }
            }
            
            for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					lowPassResult[i][j] = round(temp[i][j] / pow((weight+2),2) );
        }        

		void unsharp_masking(int b, double c){
			low_pass(b);
			int maskResult[Size][Size] = {{0}};
			
			for(int i = 0; i < Size; i++){
				for(int j = 0; j < Size; j++){
					int temp = round(ImageData[i][j]*c/(2*c-1) - lowPassResult[i][j]*(1-c)/(2*c-1));
					if(temp < 0)	temp = 0;
					maskResult[i][j] = temp;
				}
			}

			for(int i=0; i<Size; i++){
				for(int j=0; j<Size; j++){
					ImageData[i][j] = maskResult[i][j];
				}
			}
		}

		int bilinear_interpolation(float x, float y, int type){
			//type=1, for creat_circle; type=2, for create_swirl
			
			int x1 = floor(x);
			int x2 = ceil(x);
			int y1 = floor(y);
			int y2 = ceil(y);

			float dx1 = x - x1;
			float dx2 = 1 - dx1;
			float dy1 = y - y1;
			float dy2 = 1 - dy1;

			float R1 = 0, R2 =0;

			switch(type){
				case 1:
					R1 = dx2*ImageData[x1][y1] + dx1*ImageData[x2][y1];
					R2 = dx2*ImageData[x1][y2] + dx1*ImageData[x2][y2];
					break;
				case 2:
					R1 = dx2*circle[x1][y1] + dx1*circle[x2][y1];
					R2 = dx2*circle[x1][y2] + dx1*circle[x2][y2];
					break;
				default:
					cout << "wrong type" << endl;
			}
			
			int P =round(dy2*R1 + dy1*R2);
			return P;
		}

		void create_circle(){
			int radius = Size/2; 

			for(int x=0; x<Size; x++){
				for(int y=0; y<Size; y++){
					int dx = Size/2 - x;
					int dy = Size/2 - y;
					if(pow(dx,2) + pow(dy,2) <= pow(radius,2)){
						//circle[x][y] = 0;
						// int polarX = x - Size/2;
						// int polarY = y - Size/2;
						// double d = 0, theta = 0, polarXPos = 0, polarYPos = 0, xPos = 0, yPos = 0;
						// d = sqrt(pow(polarX,2) + pow(polarY,2));
						// theta = atan2(polarY, polarX) * 180 / M_PI;
						// //atan2(GC, GR) * 180 / M_PI;
						// polarXPos = d;
						// if(theta > 45 && theta < 135) 			polarYPos = (d/sin(atan2(polarY, polarX))) * cos(atan2(polarY, polarX));
						// else if(theta < -45 && theta > -135) 	polarYPos = (d/sin(atan2(polarY, polarX))) * cos(atan2(polarY, polarX));
						// else 									polarYPos = (d/cos(atan2(polarY, polarX))) * sin(atan2(polarY, polarX));
						// //result = cos ( param * PI / 180.0 );
						// //polarYPos = (d/cos(theta)) * sin(theta);
												
						// if(polarYPos < 0)		yPos = Size/2 + abs(polarYPos);
						// else if(polarYPos > 0)	yPos = Size/2 - polarYPos;
						// if(polarXPos < 0)		xPos = Size/2 + polarXPos;
						// else if(polarXPos > 0)	xPos = Size/2 - abs(polarXPos);

						//circle[x][y] = bilinear_interpolation(xPos, yPos,1);
						circle[x][y] = ImageData[x][y];
					}
				}
			}
		}

		void create_swirl(){
			int radius = Size/2; 

			for(int x=0; x<Size; x++){
				for(int y=0; y<Size; y++){
					int dx = Size/2 - x;
					int dy = Size/2 - y;
					if(pow(dx,2) + pow(dy,2) <= pow(radius,2)){

						int circleX = x -128;
						int circleY = 128 - y;

						double d = 0;
						double polarXPos = 0, polarYPos = 0, xPos = 0, yPos = 0;

						d = sqrt(pow(circleX,2) + pow(circleY,2));
												
						double angle = M_PI / 256 * d;
						polarXPos = cos(angle)*circleX - sin(angle)*circleY;
						polarYPos = sin(angle)*circleX + cos(angle)*circleY;
												
						xPos = polarXPos + Size/2;
						yPos = Size/2 - polarYPos;

						swirl[x][y] = bilinear_interpolation(xPos,yPos,2);
					}
				}
			}
			memcpy(ImageData, swirl, sizeof(unsigned char)*Size*Size);
		}

		void test_output(char out[100]){

			ofstream myfile (out);
            if (myfile.is_open())
            {
				list<double>::iterator itX;
				for(itX=checkXPos.begin(); itX!=checkXPos.end(); ++itX){
					myfile << "xPos = " << *itX << endl;
				}
                
                list<double>::iterator itY;
				for(itY=checkYPos.begin(); itX!=checkYPos.end(); ++itY){
					myfile << "yPos = " << *itY << endl;
				}
                myfile.close();
            }
            else {
                cout << "Unable to open file";
                exit(1);
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
	memset(lowPassResult, 0, sizeof(unsigned char)*Size*Size);
	memset(circle, 255, sizeof(unsigned char)*Size*Size);
	memset(swirl, 255, sizeof(unsigned char)*Size*Size);
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

	Image C, D;

	C.unsharp_masking(2, 0.65);
	C.output(argv[2]);
	D.create_circle();
	D.create_swirl();
	D.output(argv[3]);

	return 0;
}