#ifndef _MYIMG_H_
#define _MYIMG_H_
#include <memory>

class MyImg {
public:
	const int size_x, size_y;
	const float center_x, center_y;
	std::unique_ptr<unsigned char[]> inputImg;
	MyImg(std::unique_ptr<unsigned char[]> inputImg, const int& size_x, const int& size_y);
	unsigned char get(const int& i, const int& j) const;
	unsigned char get_c(const float& x, const float& y) const;
};

#endif