#include "stdafx.h"
#include "Caculator.h"

template <typename T>

Caculator<T>::Caculator()
{
}

template <typename T>
Caculator<T>::~Caculator()
{
}
template <typename T>
void Caculator<T>:: Swap( T& a,  T& b)
{
	T temp = a;
	a = b;
	b = temp;
}

// 递归算法求数组的组合
template <typename T>
void Caculator<T>::GetCombination(vector<T*> &list, T t[], int n, int m, int b[], int M)
{
	for (int i = n; i >= m; i--)
	{
		b[m - 1] = i - 1;
		if (m > 1)
		{
			GetCombination(list, t, i - 1, m - 1, b, M);
		}
		else
		{
			/*if (list == null)
			{
				list = new List<T[]>();
			}*/
			T temp[M];
			for (int j = 0; j < sizeof(b) / sizeof(b[0]); j++)
			{
				temp[j] = t[b[j]];
			}
			list.push_back(temp);
		}
	}
}


// 递归算法求排列
template <typename T>
void Caculator<T>::GetPermutation(vector<T*> &list, T t[], int startIndex, int endIndex)
{
	if (startIndex == endIndex)
	{
		/*if (list == null)
		{
			list = new List<T[]>();
		}*/
		T temp[sizeof(t) / sizeof(t[0])];
		//t.CopyTo(temp, 0);
		memcpy(temp, t, sizeof(t));
		list.push_back(temp);
	}
	else
	{
		for (int i = startIndex; i <= endIndex; i++)
		{
			Swap( t[startIndex],  t[i]);
			GetPermutation( list, t, startIndex + 1, endIndex);
			Swap( t[startIndex],  t[i]);
		}
	}

}

// 求从起始标号到结束标号的排列，其余元素不变
template <typename T>
vector<T*> Caculator<T>::GetPermutation(T t[], int startIndex, int endIndex)
{
	if (startIndex < 0 || endIndex > sizeof(t) / sizeof(t[0]) - 1)
	{
		return NULL;
	}
	vector<T*> list;
	GetPermutation( list, t, startIndex, endIndex);
	return list;
}

template <typename T>

vector<T*> Caculator<T>::GetPermutation(T t[])
{
	return GetPermutation(t, 0, t.Length - 1);
}

template <typename T>
vector<T*> Caculator<T>::GetPermutation(T t[], int n)
{
	if (n > sizeof(t)/sizeof(t[0]))
	{
		return null;
	}
	vector<T*> list
	vector<T*> c = GetCombination(t, n);
	for (int i = 0; i < c.size(); i++)
	{
		List<T*> l;
		GetPermutation( l, c[i], 0, n - 1);
	
		list.insert(list.end(),l.begin(),l.end());
	}
	return list;
}
template <typename T>
vector<T*> Caculator<T>::GetCombination(T[] t, int n)
{
	if (sizeof(t) / sizeof(t[0])< n)
	{
		return NULL;
	}
	int temp[n];
	List<T[*> list;
	GetCombination( list, t, t.Length, n, temp, n);
	return list;
}