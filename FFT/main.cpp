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
		tmp.real = obj1.real + num;
		tmp.Imaginary = obj1.Imaginary + num;
		return Complex(tmp);
	}
	friend Complex operator+(double num, Complex obj1)
	{
		return obj1 + num;
	}
	friend Complex operator-(Complex obj1, double num)
	{
		Complex tmp;
		tmp.real = obj1.real - num;
		tmp.Imaginary = obj1.Imaginary - num;
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
Complex W(int N, int k)
{
	Complex tmp;
	tmp.real = cos(2 * PI*k / N);
	tmp.Imaginary = -sin(2 * PI*k / N);
	return tmp;
}
//typedef _complex Complex;
void FFT(Complex a[], int len)
{
	for (int i = 0; i < 8; i++)
		a[i] = i;
	int counts = 0;
	int num = len - 1;
	while (num > 0)
	{
		num = num >> 1;
		counts++;
	}
	int max = pow(2, counts);
	while (len < max)
		len++;
	int *pos = new int[len];
	Complex *arr1 = new Complex[len];//���a[pos[i]]�ĸ���Ҷ�任
	Complex *arr2 = new Complex[len];//��Ż���DFT[a[pos[i]]�ĵ�ǰ�㸵��Ҷ�任
	Complex *temp = new Complex[len];
	string str;
	if (pos == NULL || arr1==NULL || arr2==NULL || temp==NULL )
		return;
	for (int i = 0; i < len; i++)
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
		pos[i] = result;
		str.clear();
	}
	//����a[pos[i]]�ĸ���Ҷ�任
	for (int i = 0; i < len; i++)
	{
		arr1[i] = a[pos[i]];
	}
	/*
	һ��Ҷ����Ϊ2^n�������������������Ϊn
	��a[pos[i]]�ĸ���Ҷ�任�У��������ö������ĵ�n�㣨��������Ϊ��1�㣩
	��ˣ��һ���Ҫѭ��for 1 to n-1����������һ�����ڵ�
	�����ٸ���Ҷ�任��ԭ�����ڣ���ֻ��Ҫ����log2(n)���������������õ�����ȫ�����ڵ㣬�Ϳ��Եõ��������x(n)�ĸ���Ҷ�任
	*/
	//���л���Ļ���arr1[i]�ĸ���Ҷ�任���õ���һ��
	for (int layer = log2(len); layer >= 1; layer--)
	{
		//ȷ��һ���
		int divide = len >> layer;
		//ȷ����ǰ���N
		int N = divide << 1;//N=divide*2
		//����������ÿ����������ӦΪgroup=pow(2,layer-1),��������㣬��Ҫ����������
		int group = pow(2,layer-1);
		//����ÿ���������еĳ��ȣ���������㣬ÿ���������еĳ���Ϊ1����len/pow(2,layer)
		int listlen = divide;
		int count1 = 0;
		int count2 = 0;
		for (int i = 0; i < group; i++)
		{
			for (int k=0; k < divide; k++)
			{
				arr2[count1] = arr1[count2] + W(N, k)*arr1[count2+divide];
				arr2[count1+len/2] = arr1[count2] - W(N, k)*arr1[count2+divide];
				count1++;
				count2 +=2;
			}
		}
		for (int i = 0; i < len; i++)
		{
			arr1[i] = arr2[i];
			cout << arr2[i] << endl;
		}
		memset(arr2, 0, 2*sizeof(double)*len);
		cout << endl;
	}
	return;
}
int main()
{
	/*Complex data[8];
	for (int i = 0; i < 8; i++)
		data[i] = i;
	FFT(data,8);*/
	/*Complex X3 = W(1, 0) * 0 + W(2, 0)*W(1, 0) * 4;
	Complex X4 = W(1, 0) * 2 + W(2, 0)*W(1, 0) * 6;
	Complex X5 = W(1, 0) * 1 + W(2, 0)*W(1, 0) * 5;
	Complex X6 = W(1, 0) * 3 + W(2, 0)*W(1, 0) * 7;
	Complex X1 = X3 + W(4, 0)*X4;
	Complex X2 = X5 + W(4, 0)*X6;
	Complex X = X1 + W(8, 0)*X2;
	cout << X << endl;*/
	int a, b;
	while (1)
	{
		cin >> a >> b;
		cout << "W(N,k): " << W(a, b) << endl;
	}
	return 0;
}