#include <vector>
#include <numeric>
#include <cstdio>
#include <cstdlib>

__int128 factorial(__int128 k) {
	if (k <= 1) {
		return 1;
	} else {
		return k * factorial(k - 1);
	}
}

std::vector<std::vector<int>> nchoosek(int n, int k) {
	std::vector<std::vector<int>> rc;
	std::vector<int> combo(k);
	std::iota(combo.begin(), combo.end(), 0);
	if (n == k) {
		rc.push_back(combo);
	} else {
		bool done = false;
		while (!done) {
			rc.push_back(combo);
			int dim = k - 1;
			while (combo[dim] == dim + n - k) {
				dim--;
				if (dim < 0) {
					done = true;
					break;
				}
			}
			if (!done) {
				combo[dim]++;
				for (int i = dim + 1; i < k; i++) {
					combo[i] = combo[i - 1] + 1;
				}
			}
		}
	}
	rc.push_back(combo);
	return rc;
}
