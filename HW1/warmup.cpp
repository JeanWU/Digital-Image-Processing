#include "stdio.h"
#include "iostream"
#include "stdlib.h"
#include "math.h"

#define Size 256
using namespace std;

unsigned char Original_img[Size][Size];

class Image{
	private:
		unsigned char ImageData[Size][Size];
	public:
		Image();
		
        void flip(){
            char tmp[Size][Size];
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    tmp[i][j] = ImageData[i][Size-1-j];
                    ImageData[i][Size-1-j] = 0;
                }
            }
            for(int i=0; i<Size; i++){
                for(int j=0; j<Size; j++){
                    ImageData[i][j] = tmp[i][j];
                }
            }
        }

		void power_trans(double power){
			double tmp[Size][Size];
            //normalize + power
			for(int i=0; i<Size; i++){
				for(int j=0; j<Size; j++){
					tmp[i][j] = pow(ImageData[i][j]/255., power);
				}
            }
            
			for(int i = 0; i < Size; i++)
				for(int j = 0; j < Size; j++)
					ImageData[i][j] = round(tmp[i][j]*255);
            
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
	for(int i = 0; i < Size; i++)
		for(int j = 0; j < Size; j++)
			ImageData[i][j] = Original_img[i][j];
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
    
	B.flip();
    B.power_trans(2.5);
	B.output(argv[2]);
    
	return 0;
}	