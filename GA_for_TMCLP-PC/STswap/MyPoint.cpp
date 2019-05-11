#include "stdafx.h"
#include "MyPoint.h"

MyPoint::MyPoint()
{
	
}

MyPoint::MyPoint(double x, double y ,int index, double STrate)
{
	this->x = x;
	this->y = y;
	this->index = index;
	this->STrate = STrate;
}


MyPoint::~MyPoint()
{
}
