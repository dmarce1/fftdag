/*
 * network.hpp
 *
 *  Created on: Mar 22, 2023
 *      Author: dmarce1
 */

#ifndef NETWORK_HPP_
#define NETWORK_HPP_

#include <vector>
#include <unordered_map>
#include <cmath>

class ann {
	int N;
	int M;
	std::vector<std::vector<std::unordered_map<int, double>>>W;
	std::vector<double> out;
	double sigmoid(double x) {
		return std::erf(x);
	}
	double rand1() {
		return (rand()/(double(RAND_MAX)+0.5));
	}
	double normal() {
		double x = rand1();
		double y = rand1();
		return sqrt(-2.0*log(x))*cos(2.0*M_PI*y);
	}
public:
	void mutate() {
		int NN = 32;
		double w0 = 0.01;
		for( int l = 0; l < M; l++) {
			for( int j = 0; j < N; j++) {
				std::vector<std::pair<int,double>> news;
				std::vector<int> kills;
				for( auto k = W[l][j].begin(); k != W[l][j].end(); k++) {
					if( rand() % NN == 0) {
						if( rand() % 4 != 0 ) {
							k->second += w0 * normal();
						} else {
							int l;
							if ( rand() % 2) {
								l = std::min(k->first+1,N-1);
							} else {
								l = std::max(0,k->first-1);
							}
							news.push_back(std::make_pair(l,w0*normal()));
						}
					} else if( rand() % (NN*10) == 0 ) {
						kills.push_back(k->first);
					}
				}
				for( auto n : news ) {
					W[l][j][n.first] = n.second;
				}
				for( auto k : kills ) {
					if(W[l][j].size()) {
						W[l][j].erase(k);
					}
				}
			}
		}
	}
	ann(int n, int m) {
		N = n;
		M = m;
		W.resize(M,std::vector<std::unordered_map<int,double>>(N));
		for( int m = 0; m < M; m++) {
			for( int j = 0; j < N; j++) {
				W[m][j][j] = 1.0;
			}
		}
		mutate();
	}
	std::vector<double> evaluate(std::vector<bool> in) {
		std::vector<double> U(N);
		std::vector<double> V(N);
		for( int n = 0; n < N; n++) {
			V[n] = in[n] ? -1.0 : 1.0;
		}
		for( int l = 0; l < M; l++) {
			for( int j = 0; j < N; j++) {
				U[j] = 0.0;
				for( auto k = W[l][j].begin(); k != W[l][j].end(); k++) {
					U[j] += V[k->first] * k->second;
				}
				U[j] = sigmoid(U[j]);
			}
			std::swap(U,V);
		}
		return std::move(V);
	}
};

#endif /* NETWORK_HPP_ */
