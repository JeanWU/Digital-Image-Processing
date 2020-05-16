#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "string.h"
#define _USE_MATH_DEFINES
#include <cmath>

 
#define Size 256
#define PI 3.14159265
#define thetaSize 181
#define rhoSize 363
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
        unsigned char gradient_map[Size][Size];
        unsigned char edge_map[Size][Size];
        unsigned char accumulate_array[thetaSize][rhoSize] = {{0}};
        unsigned char adjust_accumulate_array[thetaSize][rhoSize] = {{0}};
        unsigned char rgb_edge_map[Size][Size][3] = {{{0}}};

	public:
		Image();

        double calcuate3pointGradient(int j, int k){
            int GR = 0, GC = 0;
            if(k > 0){
                GR = ImageData[j][k] - ImageData[j][k-1];
            }
            if(j < Size - 1){
                GC = ImageData[j][k] - ImageData[j+1][k];
            }
            return sqrt(pow(GR,2) + pow(GC,2));
        }

        void computeGradientMap(){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    gradient_map[i][j] = calcuate3pointGradient(i,j);
                }
            }
        }
        
        void out_hist(char out[100]){
            int cnt[Size] = {0};
			for(int i = 0; i < Size; i++){
                for(int j = 0; j < Size; j++){
                    cnt[gradient_map[i][j]]++;
                }
            }
            //write csv file           
            ofstream myfile (out);
            if (myfile.is_open())
            {
                myfile << "value,count" << endl;
                for(int i=0; i<Size; i++){
                    myfile << i << "," << cnt[i] << endl;
                }
                myfile.close();
            }
            else {
                cout << "Unable to open file";
                exit(1);
            }
        }
        
        void generateEdgeMap(int threshold){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    if(gradient_map[i][j] >= threshold){
                        edge_map[i][j] = 255;
                    }
                }
            }
            //memcpy(ImageData, edge_map, sizeof(unsigned char)*Size*Size);
        }

        void assignAccumulateArray(){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    if(edge_map[i][j] == 255){
                        for(int theta=0; theta<=180; theta++){
                            double rho = (i-128) * cos(theta * PI / 180.0) + (j-128) * sin(theta * PI / 180.0);
                            int intRho = int(round(rho));
                            accumulate_array[theta][((rhoSize-1)/2)-intRho] += 1;
                        }
                    }
                }
            }
        }

        void adjustContrast(){
            int maxIntensity = 0;
            // find maximum number
            for(int i=0; i<thetaSize; i++){
                for(int j=0; j<rhoSize; j++){
                    if(accumulate_array[i][j] > maxIntensity){
                        maxIntensity = int(accumulate_array[i][j]);
                    }
                }
            }

            // cout << "max intensity = " << maxIntensity << endl;

            // normalize intensity 
            for(int i=0; i<thetaSize; i++){
                for(int j=0; j<rhoSize; j++){
                    int newIntensity = accumulate_array[i][j] * (255./maxIntensity);
                    adjust_accumulate_array[i][j] = newIntensity;
                }
            }
        }

        void drawSignificantLines(int LineSize){
            for(int count=0; count<LineSize; count++){
                int maxIntensity = 0;
                int rho, theta;
                // find maximum, and then assign corresponding Rho and Theta
                for(int i=0; i<thetaSize; i++){
                    // cout << "i = " << i << endl;
                    for(int j=0; j<rhoSize; j++){
                        if(adjust_accumulate_array[i][j] > maxIntensity){
                            maxIntensity = int(adjust_accumulate_array[i][j]);
                            theta = i;
                            //rho = abs(j - ((rhoSize-1)/2));
                            rho = (((rhoSize-1)/2) - j);
                        }
                    }
                }

                // cout << "theta = " << theta << "; rho = " << rho << endl;
                // cout << "jPos = " << jPos << endl;

                int x1, x2, y1, y2;
                if (theta > 180)    theta = 180;
                if (theta == 179)   theta = 180;
                if (theta == 1)     theta = 0;
                if (theta == 0 || theta == 180){
                    y1 = 0;
                    x1 = int(rho - ((y1-128)*sin(theta * PI / 180.0)) * (1./cos(theta * PI / 180.0)) ) + 128;
                    y2 = 256;
                    x2 = int(rho - ((y2-128)*sin(theta * PI / 180.0)) * (1./cos(theta * PI / 180.0)) ) + 128;
                } else if (theta == 90) {
                    x1 = 0;
                    y1 = int(rho - ((x1-128)*cos(theta * PI / 180.0)) * (1./sin(theta * PI / 180.0)) ) + 128;
                    x2 = 256;
                    y2 = int(rho - ((x2-128)*cos(theta * PI / 180.0)) * (1./sin(theta * PI / 180.0)) ) + 128;
                } else if (theta>0 && theta<90){
                    x1 = 0;
                    y1 = int(rho - ((x1-128)*cos(theta * PI / 180.0)) * (1./sin(theta * PI / 180.0)) ) + 128;
                    x2 = 256;
                    y2 = int(rho - ((x2-128)*cos(theta * PI / 180.0)) * (1./sin(theta * PI / 180.0)) ) + 128;
                }
                else {
                    x1 = 0;
                    y1 = int(rho - ((x1-128)*cos(theta * PI / 180.0)) * (1./sin(theta * PI / 180.0)) ) + 128;
                    x2 = 256;
                    y2 = int(rho - ((x2-128)*cos(theta * PI / 180.0)) * (1./sin(theta * PI / 180.0)) ) + 128;
                }

                // cout << "x1 = " << x1 << "; y1 = " << y1 << endl;
                // cout << "x2 = " << x2 << "; y2 = " << y2 << endl;

                switch(LineSize){
                    case 10:
                        drawLine_p(y1, x1, y2, x2, 'a');
                        break;
                    case 20:
                        drawLine_p(y1, x1, y2, x2, 'b');
                        break;
                    default: 
                        cout << "wrong LineSize!" << endl;
                }

                // adjust_accumulate_array[iPos][jPos] = 0;
                // clear the maximum value and its 8 neighbor
                //clear 
                if(theta == 0){
                    for(int j=-3; j<4; j++){
                        adjust_accumulate_array[0][rho+j+((rhoSize-1)/2)] = 0;
                        adjust_accumulate_array[0][abs(rho+j-((rhoSize-1)/2))] = 0;
                        adjust_accumulate_array[180][rho+j+((rhoSize-1)/2)] = 0;
                        adjust_accumulate_array[180][abs(rho+j-((rhoSize-1)/2))] = 0;
                        adjust_accumulate_array[1][rho+j+((rhoSize-1)/2)] = 0;
                        adjust_accumulate_array[1][abs(rho+j-((rhoSize-1)/2))] = 0;
                    }
                } else if(theta == 180){
                    for(int j=-3; j<4; j++){
                        adjust_accumulate_array[0][rho+j+((rhoSize-1)/2)] = 0;
                        adjust_accumulate_array[0][abs(rho+j-((rhoSize-1)/2))] = 0;
                        adjust_accumulate_array[180][rho+j+((rhoSize-1)/2)] = 0;
                        adjust_accumulate_array[180][abs(rho+j-((rhoSize-1)/2))] = 0;
                        adjust_accumulate_array[179][rho+j+((rhoSize-1)/2)] = 0;
                        adjust_accumulate_array[179][abs(rho+j-((rhoSize-1)/2))] = 0;
                    }
                } else {
                    for(int i=-3; i<4; i++){
                        for(int j=-3; j<4; j++){
                            adjust_accumulate_array[theta+i][rho+j+((rhoSize-1)/2)] = 0;
                            adjust_accumulate_array[theta+i][abs(rho+j-((rhoSize-1)/2))] = 0;
                        }
                    }
                }
            }
            
        }

        void copyRgbEdgeMap(){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    for(int k=0; k<3; k++){
                        if(edge_map[i][j] == 255){
                            rgb_edge_map[i][j][k] = 255;
                        }
                    }
                }
            }
        }

        void clearRbgEdgeMap(){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    for(int k=0; k<3; k++){
                        rgb_edge_map[i][j][k] = 0;
                    }
                }
            }
        }

        void clearAdjustAccumulateArray(){
            //adjust_accumulate_array[thetaSize][rhoSize]
            for(int i=0; i<thetaSize; i++){
                for(int j=0; j<rhoSize; j++){
                    adjust_accumulate_array[i][j] = 0;
                }
            }
        }

        static inline float fastAtan2f(float dy, float dx){
            // 快速atan運算
            static const float atan2_p1 = 0.9997878412794807f*(float)(180/M_PI);
            static const float atan2_p3 = -0.3258083974640975f*(float)(180/M_PI);
            static const float atan2_p5 = 0.1555786518463281f*(float)(180/M_PI);
            static const float atan2_p7 = -0.04432655554792128f*(float)(180/M_PI);
            static const float atan2_DBL_EPSILON = 2.2204460492503131e-016;

            float ax = std::abs(dx), ay = std::abs(dy);
            float a, c, c2;
            if (ax >= ay) {
                c = ay/(ax + static_cast<float>(atan2_DBL_EPSILON));
                c2 = c*c;
                a = (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
            } else {
                c = ax/(ay + static_cast<float>(atan2_DBL_EPSILON));
                c2 = c*c;
                a = 90.f - (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
            }
            if (dx < 0)
                a = 180.f - a;
            if (dy < 0)
                a = 360.f - a;
            return a;
        }

        void drawLine_p(int y, int x, int y2, int x2, char type) {
            // 兩點之間的距離差
            float dx = static_cast<float>(x2-x);
            float dy = static_cast<float>(y2-y);
            // 以Y軸為主
            float sita=fastAtan2f(dy, dx);
            if ((sita>45 && sita<135) || (sita>225 && sita<315)) {
                float slopeY = dx/dy; // 斜率
                for (int i = 0; i < abs(dy); i++) {
                    int iFix = dy>0? i:-i;
                    int currPos = static_cast<int>(iFix*slopeY+.5f + x);

                    int distX = currPos;
                    int distY = y+iFix;

                    if (distX<0 or distX>=Size or distY<0 or distY>=Size) {
                        return;
                    }
                    //img.raw_img[distY*img.width + distX] = static_cast<unsigned char>(val);
                    switch(type){
                        // draw red line for D1
                        case 'a':
                            rgb_edge_map[distX][distY][0] = 255;
                            rgb_edge_map[distX][distY][1] = 0;
                            rgb_edge_map[distX][distY][2] = 0;
                            break;
                        // draw green line for D2
                        case 'b':
                            rgb_edge_map[distX][distY][0] = 0;
                            rgb_edge_map[distX][distY][1] = 255;
                            rgb_edge_map[distX][distY][2] = 0;
                            break;
                        default: 
                            cout << "wrong type format for drawing line!" << endl;
                    }
                    
                }
            } 
            // 以X軸為主
            else {
                float slopeX = dy/dx; // 斜率
                for (int i = 0; i < abs(dx); i++) {
                    int iFix = dx>0? i:-i;
                    int currPos = static_cast<int>(iFix*slopeX+.5 + y);

                    int distX = x+iFix;
                    int distY = currPos;

                    if (distX<0 or distX>=Size or distY<0 or distY>=Size) {
                        return;
                    }
                    
                    switch(type){
                        // draw red line for D1
                        case 'a':
                            rgb_edge_map[distX][distY][0] = 255;
                            rgb_edge_map[distX][distY][1] = 0;
                            rgb_edge_map[distX][distY][2] = 0;
                            break;
                        // draw green line for D2
                        case 'b':
                            rgb_edge_map[distX][distY][0] = 0;
                            rgb_edge_map[distX][distY][1] = 255;
                            rgb_edge_map[distX][distY][2] = 0;
                            break;
                        default: 
                            cout << "wrong type format for drawing line!" << endl;
                    }
                }
            }
        }


		void output(char out[100], char type){
			FILE *fout;
			if (!(fout=fopen(out,"wb"))){
				cout<<"Cannot open file!"<<endl;
				exit(1);
			}
            switch(type){
                //type == 'a', output edge_map
                //type == 'b', output accumulate_array
                //type == 'c', output adjust_accumulate_array
                case 'a':
                    fwrite(edge_map, sizeof(unsigned char), Size*Size, fout);
                    break;
                case 'b':
                    fwrite(accumulate_array, sizeof(unsigned char), thetaSize*rhoSize, fout);
                    break;
                case 'c':
                    fwrite(adjust_accumulate_array, sizeof(unsigned char), thetaSize*rhoSize, fout);
                    break;
                case 'd':
                    fwrite(rgb_edge_map, sizeof(unsigned char), Size*Size*3, fout);
                    break;
                default: 
                    cout << "please select correct output image type!" << endl;
            }
			
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
    
    Image sample1;
    sample1.computeGradientMap();
    //sample1.out_hist(argv[2]);
    sample1.generateEdgeMap(40);
    char E[] = "E.raw";
    sample1.output(E, 'a'); //-> E
    
    sample1.assignAccumulateArray();
    char H1[] = "H1.raw";
    sample1.output(H1, 'b'); //-> H1

    sample1.adjustContrast();
    char H2[] = "H2.raw";
    sample1.output(H2, 'c'); // -> H2

    sample1.copyRgbEdgeMap();
    sample1.drawSignificantLines(10);
    char D1[] = "D1.raw";
    sample1.output(D1, 'd');

    sample1.clearAdjustAccumulateArray();
    sample1.adjustContrast();
    sample1.clearRbgEdgeMap();
    sample1.copyRgbEdgeMap();
    sample1.drawSignificantLines(20);
    char D2[] = "D2.raw";
    sample1.output(D2, 'd');




	return 0;
}
