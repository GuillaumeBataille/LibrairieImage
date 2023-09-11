#ifndef COLOR_H
#define COLOR_H

#include <iostream>
#include "image_ppm.h"

class Color{

	public:
		int r,g,b;

        Color(){}

		Color(int r, int g, int b)
		{
			this->r = r;
			this->g = g;
			this->b = b;
		}

        Color(OCTET r, OCTET g, OCTET b){
            this->r = r;
            this->g = g;
            this->b = b;
        }

        Color& operator = (Color &c){
            r = c.r;
            g = c.g;
            b = c.b;
            return *this;
        }
        Color& operator = (const Color &c){
            r = c.r;
            g = c.g;
            b = c.b;
            return *this;
        }

        Color& operator + (Color &c){
            r = r + c.r;
            g = g + c.g;
            b = b + c.b;
            return *this;
        }

        Color& operator / (int i)
        {
            r /= i;
            g /= i;
            b /= i;
            return *this;
        }

        bool operator == (Color &c)
        {
            if( r == c.r &&
                g == c.g &&
                b == c.b )
                return true;
            return false;
        }
			
		float dist(Color c)
		{
			return sqrt(pow(r - c.r,2) + pow(g - c.g,2) + pow(b - c.b,2));
		}

        std::string toString()
        {
            return std::to_string(r) + " " + std::to_string(g) + " " + std::to_string(b);
        }
};

class Color_float{

public:
    float r,g,b;

    Color_float(float r, float g, float b)
    {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    Color_float(){
        r = 0;
        g = 0;
        b = 0;
    }

    Color_float& operator = (Color_float &c){
        r = c.r;
        g = c.g;
        b = c.b;
        return *this;
    }
    Color_float& operator = (Color &c){
        r = c.r;
        g = c.g;
        b = c.b;
        return *this;
    }

    Color_float& operator + (Color_float &c){
        r = r + c.r;
        g = g + c.g;
        b = b + c.b;
        return *this;
    }
    Color_float& operator + (Color &c){
        r = r + c.r;
        g = g + c.g;
        b = b + c.b;
        return *this;
    }

    Color_float& operator / (int i)
    {
        r /= static_cast<float>(i);
        g /= static_cast<float>(i);
        b /= static_cast<float>(i);
        return *this;
    }

    Color convertToInt() const{
        Color c((int)r,(int)g,(int)b);
        return c;
    }

    float dist(Color_float c)
    {
        return sqrt(pow(r - c.r,2) + pow(g - c.g,2) + pow(b - c.b,2));
    }
    float dist(Color c)
    {
        return sqrt(pow(r - c.r,2) + pow(g - c.g,2) + pow(b - c.b,2));
    }

    std::string toString()
    {
        return std::to_string(r) + " " + std::to_string(g) + " " + std::to_string(b);
    }
};

#endif
