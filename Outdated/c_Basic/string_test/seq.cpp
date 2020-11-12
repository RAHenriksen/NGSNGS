#include <iostream>
#include <string>
#include <vector>

vector<int> findLocation(const string sample, char findIt){
	vector<int> characterLocations;
	    for(int i =0; i < sample.size(); i++)
			if(sample[i] == findIt)
				characterLocations.push_back(i);
				return characterLocations;
}

int main() {
	std::string DNA = "ATGGTACCTAG";
	std::cout << "DNA WAS LOADED\n" << DNA << std::endl;
	std::cout << DNA.length() << std::endl;
	// std::string DNA2(DNA.begin(), std::find(DNA.begin(), DNA.end(), 'T'));
    // std::cout << DNA2;
	std::cout << findLocation("GGTATAC",'T');
	return 0;
}
