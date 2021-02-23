#include "MyImg.h"

MyImg::MyImg(std::unique_ptr<unsigned char[]> i_inputImg, const int& i_size_x, const int& i_size_y) 
	: size_x(i_size_x), size_y(i_size_y), inputImg(std::move(i_inputImg)), center_x(i_size_x * 0.5f), center_y(i_size_y * 0.5f)
{};

unsigned char MyImg::get(const int& i, const int& j) const{
	if (i < 0 || i >= size_x || j < 0 || j >= size_y) {
		return 0;
	}
	return inputImg[i + j * size_x];
};

unsigned char MyImg::get_c(const float& x, const float& y) const {
	if (x < -center_x || x >= center_x || y < -center_y || y >= center_y) {
		return 0;
	}
	return inputImg[floor((float)x + center_x) + floor((float)y + center_y) * size_x];
};