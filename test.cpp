#include <vector>
#include <iostream>

using namespace std;

int main(){
	vector <int> v;

	cout << v.capacity() << " " << v.size() << endl;

	v.push_back(22);
	cout << v.capacity() << " " << v.size() << endl;

	v.reserve(100);
	cout << v.capacity() << " " << v.size() << endl;
}
