#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#define PI 3.141592653589793238426
using namespace std;
class Complex
{
public:
	double real;
	double Imaginary;
	Complex()
	{
		this->real = 0;
		this->Imaginary = 0;
	}
	Complex(double x, double y)
	{
		this->real = x;
		this->Imaginary = y;
	}
	Complex(const Complex &obj)
	{
		this->real = obj.real;
		this->Imaginary = obj.Imaginary;
	}
	Complex operator=(double x)
	{
		real = x;
		Imaginary = 0.0;
		return *this;
	}
	friend Complex operator+(Complex obj1, double num)
	{
		Complex tmp;
		tmp.real = obj1.real+num;
		tmp.Imaginary = obj1.Imaginary+num;
		return Complex(tmp);
	}
	friend Complex operator+(double num, Complex obj1)
	{
		return obj1 + num;
	}
	friend Complex operator-(Complex obj1, double num)
	{
		Complex tmp;
		tmp.real = obj1.real-num;
		tmp.Imaginary = obj1.Imaginary-num;
		return Complex(tmp);
	}
	friend Complex operator-(double num, Complex obj1)
	{
		return obj1 + (-num);
	}
	friend Complex operator+(Complex obj1, Complex obj2)
	{
		Complex tmp;
		tmp.real = obj1.real + obj2.real;
		tmp.Imaginary = obj1.Imaginary + obj2.Imaginary;
		return Complex(tmp);
	}
	friend Complex operator-(Complex obj1, Complex obj2)
	{
		Complex tmp;
		tmp.real = obj1.real - obj2.real;
		tmp.Imaginary = obj1.Imaginary - obj2.Imaginary;
		return Complex(tmp);
	}
	friend Complex operator*(Complex obj1, double num)
	{
		Complex tmp;
		tmp.real = obj1.real*num;
		tmp.Imaginary = obj1.Imaginary*num;
		return Complex(tmp);
	}
	friend Complex operator*(double num, Complex obj1)
	{
		return obj1 * num;
	}
	friend Complex operator*(Complex obj1, Complex obj2)
	{
		Complex tmp;
		tmp.real = obj1.real*obj2.real - obj1.Imaginary*obj2.Imaginary;
		tmp.Imaginary = obj1.real*obj2.Imaginary + obj1.Imaginary*obj2.real;
		return Complex(tmp);
	}
	friend Complex operator/(Complex obj1, double num)
	{
		_STL_VERIFY((num != 0), "Division cannot be zero!");
		Complex tmp;
		tmp.real = obj1.real / num;
		tmp.Imaginary = obj1.Imaginary / num;
		return Complex(tmp);
	}
	friend Complex operator/(Complex obj1, Complex obj2)
	{
		_STL_VERIFY((obj2.real != obj2.Imaginary), "Division cannot be zero!");
		Complex tmp;
		Complex tmp2 = obj2;
		tmp2.Imaginary = -obj2.Imaginary;
		tmp = (obj1 * tmp2) / (obj2.real*obj2.real + obj2.Imaginary + obj2.Imaginary);
		return Complex(tmp);
	}
	friend ostream &operator<<(ostream &os, const Complex obj)
	{
		double epsilon = 1e-4;
		if (abs(obj.real - 0.0) < epsilon)
		{
			if (abs(obj.Imaginary - -1.0) < epsilon)
				os << "-i";
			else if (abs(obj.Imaginary - 1.0) < epsilon)
				os << "i";
			else if (abs(obj.Imaginary - 0.0) < epsilon)
				os << fixed << setprecision(3) << abs(obj.real);
			else
				os << fixed << setprecision(3) << obj.Imaginary << "i";
		}
		else
		{
			if (abs(obj.Imaginary - -1.0) < epsilon)
				os << fixed << setprecision(3) << obj.real << " - i";
			else if (abs(obj.Imaginary - 1.0) < epsilon)
				os << fixed << setprecision(3) << obj.real << " + i";
			else if (abs(obj.Imaginary - 0.0) < epsilon)
				os << fixed << setprecision(3) << obj.real;
			else if (obj.Imaginary < -epsilon)
				os << fixed << setprecision(3) << obj.real << " - " << setprecision(3) << abs(obj.Imaginary) << "i";
			else
				os << fixed << setprecision(3) << obj.real << " + " << setprecision(3) << obj.Imaginary << "i";
		}
		return os;
	}
	friend istream &operator>>(istream &is, Complex obj)
	{
		is >> obj.real >> obj.Imaginary;
		return is;
	}
};
double sqrt(Complex obj)
{
	return std::sqrt(obj.real*obj.real + obj.Imaginary*obj.Imaginary);
}
Complex W(int N,int k)
{
	Complex tmp;
	tmp.real = cos(2 * PI*k/N);
	tmp.Imaginary = -sin(2 * PI*k / N);
	return tmp;
}
//typedef _complex Complex;
void FFT(Complex a[],int len)
{
	int *pos=new int[len];
	string str;
	if (pos == NULL)
		return;
	int counts = 0;
	int num = len;
	while (num > 0)
	{
		num /= 2;
		counts++;
	}
	for (int i = 0; i <= len; i++)
	{
		int count = i;
		int result = 0;
		if (count == 0)
			str.append("0");
		while (count > 0)
		{
			int tmp = count % 2;
			if (tmp == 0)
				str.append("0");
			else
				str.append("1");
			count = count >> 1;
		}
		while (str.length() < counts)
			str.append("0");
		reverse(str.begin(), str.end());
		for (int i = 0; i < str.length(); i++)
			result += (str[i] - 48)*pow(2, i);
		a[i] = result;
		str.clear();
	}

	return;
}
int main()
{
	Complex data[7];
	for (int i = 0; i < 7; i++)
		data[i] = i;
	//FFT(data,7);
	Complex X3 = W(1, 0) * 0 + W(2, 0)*W(1, 0) * 4;
	Complex X4 = W(1, 0) * 2 + W(2, 0)*W(1, 0) * 6;
	Complex X5 = W(1, 0) * 1 + W(2, 0)*W(1, 0) * 5;
	Complex X6 = W(1, 0) * 3 + W(2, 0)*W(1, 0) * 7;
	Complex X1 = X3 + W(4, 0)*X4;
	Complex X2 = X5 + W(4, 0)*X6;
	Complex X = X1 + W(8, 0)*X2;
	cout<<X<<endl;
	return 0;
}
