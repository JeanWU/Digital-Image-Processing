#include "stdio.h"
#include "iostream"
#include "fstream"
#include "stdlib.h"
#include "string.h"
 
#define Size 256
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
        unsigned char Gbefore[Size][Size] = {{0}};
        unsigned char Gafter[Size][Size] = {{0}};
        unsigned char label[Size][Size] = {{0}};
        unsigned char rgbImageData[Size][Size][3] = {{{0}}};
        int find = 1;

	public:
		Image();

        
        void boundary_extraction(){
            unsigned char G[Size][Size] = {{0}};

            for(int x=2; x<Size-2; x++){
                for(int y=2; y<Size-2; y++){
                    if( erosion(x,y,1) == 1)  G[x][y] = 1;
                }
            }

            for(int x=2; x<Size-2; x++){
                for(int y=2; y<Size-2; y++){
                    if( G[x][y] == 1)  ImageData[x][y] = 0;
                }
            }
        }

        int erosion(int x, int y, int maskType){
            //int mask[5][5] = {{0,0,1,0,0},{0,1,1,1,0},{1,1,1,1,1},{0,1,1,1,0},{0,0,1,0,0}};
            int mask1[3][3] = {{0,1,0},{1,1,1},{0,1,0}};
            int mask2[3][3] = {{1,1,1},{1,1,1},{1,1,1}};
            int maskSize = 3;
            int offset = (maskSize - 1) / 2;
            for(int i=0; i<maskSize; i++){
                for(int j=0; j<maskSize; j++){
                    switch(maskType){
                        case 1: 
                            if(mask1[i][j] == 1 && ImageData[x+i-offset][y+j-offset] == 0)     return 0;
                            break;
                        case 2: 
                            if(mask2[i][j] == 1 && ImageData[x+i-offset][y+j-offset] == 0)     return 0;
                            break;
                        default:    
                            cout << "wrong mask type!" << endl;
                    }
                    
                }
            }
            return 1;
        }

        void connect_component_labeling(){           
            for(int x=1; x<Size-1; x++){
                for(int y=1; y<Size-1; y++){
                    if(ImageData[x][y] == 255 && label[x][y] == 0){
                        Gbefore[x][y] = 255;
                        connect();                        
                    }
                }
            }
            drawlabel();
        }

        void drawlabel(){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    if(label[i][j] == 1){
                        rgbImageData[i][j][0] = 255;
                        rgbImageData[i][j][1] = 0;
                        rgbImageData[i][j][2] = 0;
                    }
                    if(label[i][j] == 2){
                        rgbImageData[i][j][0] = 0;
                        rgbImageData[i][j][1] = 255;
                        rgbImageData[i][j][2] = 0;
                    }
                    if(label[i][j] == 3){
                        rgbImageData[i][j][0] = 0;
                        rgbImageData[i][j][1] = 0;
                        rgbImageData[i][j][2] = 255;
                    }
                    if(label[i][j] == 4){
                        rgbImageData[i][j][0] = 255;
                        rgbImageData[i][j][1] = 255;
                        rgbImageData[i][j][2] = 0;
                    }
                    if(label[i][j] == 5){
                        rgbImageData[i][j][0] = 0;
                        rgbImageData[i][j][1] = 255;
                        rgbImageData[i][j][2] = 255;
                    }
                    if(label[i][j] == 6){
                        rgbImageData[i][j][0] = 255;
                        rgbImageData[i][j][1] = 0;
                        rgbImageData[i][j][2] = 255;
                    }
                    if(label[i][j] == 7){
                        rgbImageData[i][j][0] = 255;
                        rgbImageData[i][j][1] = 255;
                        rgbImageData[i][j][2] = 255;
                    }
                    if(label[i][j] == 8){
                        rgbImageData[i][j][0] = 128;
                        rgbImageData[i][j][1] = 128;
                        rgbImageData[i][j][2] = 128;
                    }
                }
            }

        }

        void connect(){
            for(int x=1; x<Size-1; x++){
                for(int y=1; y<Size-1; y++){
                    if (Gbefore[x][y] == 255){
                        dilation(x, y);
                    }
                }
            }
            if (compareBeforeAfter() == false){
                memcpy(Gbefore, Gafter, sizeof(unsigned char)*Size*Size);
                connect();
            }
            else if (compareBeforeAfter() == true){
                //memcpy(ImageData, Gafter, sizeof(unsigned char)*Size*Size);
                updatelabel();
                updateG();
            }
        }

        void updatelabel(){
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    if(Gafter[i][j] == 255) label[i][j] = find;
                }
            }
            find ++;
        }

        void updateG(){
            memset(Gbefore, 0, Size*Size*sizeof(unsigned char));
            memset(Gafter, 0, Size*Size*sizeof(unsigned char));
        }

        void dilation(int x, int y){
            int H[3][3] = {{1,1,1}, {1,1,1}, {1,1,1}};
            int Hsize = 3;
            int Hoffset = (Hsize - 1) / 2;
            for(int i=0; i<Hsize; i++){
                for(int j=0; j<Hsize; j++){
                    if(H[i][j] == 1 && ImageData[x+i-Hoffset][y+j-Hoffset] == 255){
                        Gafter[x+i-Hoffset][y+j-Hoffset] = 255;
                    }
                }
            }
        }

        int dilation(int x, int y, unsigned char Image[][256]){
            int H[3][3] = {{1,1,1}, {1,1,1}, {1,1,1}};
            int Hsize = 3;
            int Hoffset = (Hsize - 1) / 2;
            for(int i=0; i<Hsize; i++){
                for(int j=0; j<Hsize; j++){
                    if(H[i][j] == 1 && Image[x+i-Hoffset][y+j-Hoffset] == 255){
                        return 1;
                    }
                }
            }
            return 0;
        }

        bool compareBeforeAfter(){
            for(int x=1; x<Size-1; x++){
                for(int y=1; y<Size-1; y++){
                    if (Gbefore[x][y] != Gafter[x][y])
                        return false;
                }
            }
            return true;
        }

        void thinning(){
            int S = 0, N = 0;
            unsigned char mark[Size][Size] = {{0}};
            bool flag = true;

            while(flag){
                for(int step=1; step<3; step++){
                    for(int x=1; x<Size-1; x++){
                        for(int y=1; y<Size-1; y++){
                            S = count_transition_neighbor(x,y,1);
                            N = count_transition_neighbor(x,y,2);
                            if( (ImageData[x][y] == 255) && (N>=2) && (N<=6) && (S==1) && (checkP(x,y,step)==true) ){
                                mark[x][y] = 1;
                            }
                        }
                    }
                    int countMark = 0;
                    for(int x=1; x<Size-1; x++){
                        for(int y=1; y<Size-1; y++){
                            if (mark[x][y] == 1){
                                ImageData[x][y] = 0;
                                countMark ++;
                            }
                        }
                    }
                    if(countMark == 0){
                        flag = false;
                    }
                    memset(mark, 0, Size*Size*sizeof(unsigned char));
                }
            }
        }

        bool checkP(int x, int y, int type){
            //if type = 1, check P4=0 or P6=0 or (P2=P8=0)
            //if type = 2, check P2=0 or P8=0 or (P4=P6=0)
            int sequence[8][2] = {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}};
            signed char P[8] = {0};
            for(int k=0; k<8; k++){
                int i = sequence[k][0];
                int j = sequence[k][1];
                if(ImageData[x+i][y+j] == 255)  P[k] = 1;
                if(ImageData[x+i][y+j] == 0)  P[k] = 0;
            }

            switch(type){
                case 1:
                    if((P[2] == 0) || (P[4] == 0) || (P[0]==P[6]==0)){
                        return true;
                    }
                    else{
                        return false;
                    }
                    break;
                case 2: 
                    if((P[0] == 0) || (P[6] == 0) || (P[2]==P[4]==0)){
                        return true;
                    }
                    else{
                        return false;
                    }
                    break;
                default:
                    cout << "wrong type!" << endl;
                    return false;
            }
        }

        int count_transition_neighbor(int x, int y, int type){
            int sequence[8][2] = {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}};
            signed char P[8] = {0};
            int countTransition = 0;
            int countNeighbor = 0;
            for(int k=0; k<8; k++){
                int i = sequence[k][0];
                int j = sequence[k][1];
                if(ImageData[x+i][y+j] == 255)  P[k] = 1;
                if(ImageData[x+i][y+j] == 0)  P[k] = 0;
            }

            switch(type){
                case 1: //count transition
                    for(int k=0; k<8; k++){
                        if(P[k]==0 && P[k+1]==1){
                            countTransition ++;
                        }
                    }
                    return countTransition;
                    break;
                case 2: //count neighbor
                    for(int k=0; k<8; k++){
                        if(P[k] == 1){
                            countNeighbor ++;
                        }
                    }
                    return countNeighbor;
                    break;
                default:
                    cout << "wrong type!" << endl;
                    return 0;
            }
        }

        void skeleton(){
            unsigned char temp[Size][Size] = {{0}};
            unsigned char E[7][Size][Size] = {{{0}}};
            unsigned char O[7][Size][Size] = {{{0}}};
            unsigned char S[7][Size][Size] = {{{0}}};
            
            //erosion => E0, E1, E2, E3, E4, E5, E6
            for(int k=0; k<7; k++){
                for(int x=1; x<Size-1; x++){
                    for(int y=1; y<Size-1; y++){
                        if( erosion(x,y,2) == 1)  temp[x][y] = 1;
                    }
                }
                for(int x=1; x<Size-1; x++){
                    for(int y=1; y<Size-1; y++){
                        if( temp[x][y] == 0)  ImageData[x][y] = 0;
                    }
                }
                memcpy(E[k], ImageData, sizeof(unsigned char)*Size*Size);
                memset(temp, 0, sizeof(unsigned char)*Size*Size);
            }
            //opening => O0, O1, ..., O6
            //O6 should be null by default
            for(int k=0; k<6; k++){
                for(int x=1; x<Size-1; x++){
                    for(int y=1; y<Size-1; y++){
                        if( dilation(x, y, E[k+1]) == 1)  temp[x][y] = 1;
                    }
                }
                for(int x=1; x<Size-1; x++){
                    for(int y=1; y<Size-1; y++){
                        if( temp[x][y] == 1)  O[k][x][y] = 255;
                    }
                }
                memset(temp, 0, sizeof(unsigned char)*Size*Size);
            }
            //Skeleton => S0, S1, ..., S6
            //S0 = E0 - O0
            for(int k=0; k<7; k++){
                for(int x=1; x<Size-1; x++){
                    for(int y=1; y<Size-1; y++){
                        S[k][x][y] = E[k][x][y] - O[k][x][y];
                    }
                }
            }
            //S(A) = S0+S1+...+S6; store the result in ImageData
            for(int k=0; k<7; k++){
                for(int x=1; x<Size-1; x++){
                    for(int y=1; y<Size-1; y++){
                        if(S[k][x][y] == 255) ImageData[x][y] = 255;
                    }
                }
            }           
        }

        bool is_null(){
            int count255 = 0;
            for(int x=1; x<Size-1; x++){
                for(int y=1; y<Size-1; y++){
                    if(ImageData[x][y] == 255) count255++;
                }
            }
            if(count255 == 0)   return true;
            else                return false;
        }
        		
		void output(char out[100], int type){
			FILE *fout;
			if (!(fout=fopen(out,"wb"))){
				cout<<"Cannot open file!"<<endl;
				exit(1);
			}
            switch(type){
                case 1:
                    fwrite(ImageData, sizeof(unsigned char), Size*Size, fout);
                    break;
                case 2:
                    fwrite(rgbImageData, sizeof(unsigned char), Size*Size*3, fout);
                    break;
                default:
                    cout << "wrong type!" << endl;
            }           
			fclose(fout);
		}
        
};
Image::Image(){
	memcpy(ImageData, Original_img, sizeof(unsigned char)*Size*Size);
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
    
    Image B;
    B.boundary_extraction();
    B.output(argv[2], 1);
    Image C;
    C.connect_component_labeling();
    C.output(argv[3], 2);
    Image D1;
    D1.thinning();
    D1.output(argv[4], 1);
    Image D2;
    D2.skeleton();
    D2.output(argv[5], 1);

	return 0;
}