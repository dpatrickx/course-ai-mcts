#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include "Point.h"
#include "Judge.h"
#include "Strategy.h"
using namespace std;

#define NODENUM 10000000
#define epsilon 1e-6
#define C 2.0
#define TIMELIMIT 4.5

class MCNode
{
public:
	int x,y;
	int lChild,rChild,father;
	bool user,isLeaf;				// user: false-->oppoent, true-->me
	int winRound, totRound;

	MCNode() : x(0),y(0),lChild(-1),rChild(-1),father(-1),user(false),isLeaf(true),winRound(0),totRound(0)
	{}
};

MCNode nodes[NODENUM];
int rank;
int sel;
int fNode = 0;
int* tempTop1;
int* tempTop2;				// used for modify function
int** tempBoard1;
int** tempBoard2;			// used for modify function
int nn = 0,mm = 0,banX = -1,banY = -1;
clock_t start,finish;
double totaltime = 0;

extern "C" __declspec(dllexport) Point* getPoint(const int M, const int N, const int* top, const int* _board, 
	const int lastX, const int lastY, const int noX, const int noY){

	int x = -1, y = -1;//最终将你的落子点存到x,y中
	int** board = new int*[M];
	for(int i = 0; i < M; i++){
		board[i] = new int[N];
		for(int j = 0; j < N; j++){
			board[i][j] = _board[i * N + j];
		}
	}
// start of my part
	srand(unsigned(time(0)));
	start = clock();
	totaltime = 0;
	// time inited
	nn = N;
	mm = M;
	banX = noX;
	banY = noY;
	tempTop1 = new int[N];
	tempTop2 = new int[N];
	tempBoard1 = new int* [M];
	tempBoard2 = new int* [M];
	for (int i = 0;i < M;i++){
		tempBoard1[i] = new int[N];
		tempBoard2[i] = new int[N];
		for (int j = 0;j < N;j++){
			tempBoard1[i][j] = board[i][j];
			tempBoard2[i][j] = board[i][j];
		}
	}
	for (int i = 0;i < N;i++){
		tempTop1[i] = top[i];
		tempTop2[i] = top[i];
	}
	nodes[0].isLeaf = true;
	nodes[0].user = false;
	fNode = 0;
	rank = 1;
	nodes[0].totRound = 0;
	// parameters inited
	while(totaltime < TIMELIMIT){
		if (nodes[fNode].isLeaf){
			// expands begins
			if (tempTop1[banY]-1 == banX)
				tempTop1[banY]--;
			// ensure every node in nodes array is available
			nodes[fNode].lChild = rank;
			for (int i = 0;i < N;i++){
				if (tempTop1[i] > 0){
					nodes[rank].father = fNode;
					nodes[rank].user = !nodes[fNode].user;
					nodes[rank].x = tempTop1[i]-1;
					nodes[rank].y = i;
					nodes[rank].totRound = 0;
					nodes[rank].winRound = 0;			
					nodes[rank].isLeaf = true;
					rank++;
				}
			}
			nodes[fNode].rChild = rank;
			if (rank == nodes[fNode].lChild){
				for (int j = fNode;j != 0;j = nodes[j].father){
				nodes[j].totRound++;
				if (!nodes[j].user)
					nodes[j].winRound++;
				}
				nodes[0].totRound++;

				for (int i = 0;i < M;i++)
					for (int j = 0;j < N;j++)
					tempBoard1[i][j] = board[i][j];
				for (int i = 0;i < N;i++)
					tempTop1[i] = top[i];
				fNode = 0;

				finish = clock();
				totaltime = ((double)finish-start)/CLOCKS_PER_SEC;				
				continue;
			}
			nodes[fNode].isLeaf = false;
			// expand ends

			// simulate begins	
			for (int i = nodes[fNode].lChild;i < nodes[fNode].rChild;i++){
				bool resultSim;
				for (int ii = 0;ii < N;ii++)
					tempTop2[ii] = tempTop1[ii];
				for (int ii = 0;ii < M;ii++)
					for (int j = 0;j < N;j++)
						tempBoard2[ii][j] = tempBoard1[ii][j];
				bool sideSim = nodes[i].user;
				int mWinSim = 0,uWinSim = 0;
				int ySim = nodes[i].y;
				int tempXSim = nodes[i].x;

				if (sideSim)
					tempBoard2[tempXSim][ySim] = 2;
				else
					tempBoard2[tempXSim][ySim] = 1;
				tempTop2[ySim] = tempXSim;
				while(1){
					// judge is the game is over
					if (sideSim)
						mWinSim = machineWin(tempXSim,ySim,mm,nn,tempBoard2);
					else
						uWinSim = userWin(tempXSim,ySim,mm,nn,tempBoard2);

					if (mWinSim)				{resultSim = 1;break;}
					if (uWinSim)				{resultSim = 0;break;}
					if (isTie(nn,tempTop2))		{resultSim = 0;break;}
					// ensure every element of tempTop2 is available
					if (tempTop2[banY]-1 == banX)
						tempTop2[banY]--;
					// select a positon available
					ySim = rand()%nn;
					while(tempTop2[ySim] == 0)
						ySim = rand()%nn;
					// update top and board
					sideSim = !sideSim;
					tempTop2[ySim]--;
					tempXSim = tempTop2[ySim];
					if (sideSim)
						tempBoard2[tempXSim][ySim] = 2;			// myself(machine)
					else
						tempBoard2[tempXSim][ySim] = 1;			// oppoent(user)
				}
				// simulate(i) ends
				// update totRound and winRound
				for (int j = i;j != 0;j = nodes[j].father){
					nodes[j].totRound++;
					if (resultSim && nodes[j].user)
						nodes[j].winRound++;
					if (!resultSim && !nodes[j].user)
						nodes[j].winRound++;
				}
				nodes[0].totRound++;
			} // end for
		} // end ifLeaf
		//******************************************
		//select() begins
		sel = nodes[fNode].lChild;
		double mSel =  ((double)nodes[sel].winRound/((double)nodes[sel].totRound+epsilon) +     
		C*sqrt((double)log((double)nodes[fNode].totRound+1)/((double)nodes[sel].totRound+epsilon)));
		for (int i = nodes[fNode].lChild+1;i < nodes[fNode].rChild;i++){
			double tempValue = ((double)nodes[i].winRound/((double)nodes[i].totRound+epsilon) +     
					C*sqrt((double)log((double)nodes[fNode].totRound+1)/((double)nodes[i].totRound+epsilon)));
			if (tempValue > mSel){
				mSel = tempValue;
				sel = i;
			}
		}
		int tempY = nodes[sel].y;
		tempTop1[tempY] = nodes[sel].x;

		if (nodes[sel].user)
			tempBoard1[tempTop1[tempY]][tempY] = 2;
		else
			tempBoard1[tempTop1[tempY]][tempY] = 1;
		
		if (tempTop1[banY]-1 == banX)
			tempTop1[banY]--;
		//select() ends
		//******************************************		

		// judge if game is over
		bool gameOver = false;
		if (nodes[sel].user){
			if (machineWin(nodes[sel].x,nodes[sel].y,M,N,tempBoard1)){
				for (int j = sel;j != 0;j = nodes[j].father){
					nodes[j].totRound++;
					if (nodes[j].user)
						nodes[j].winRound++;
				}
				nodes[0].totRound++;
				gameOver = true;
			}
		}

		if (!nodes[sel].user){
			if (userWin(nodes[sel].x,nodes[sel].y,M,N,tempBoard1)){
				for (int j = sel;j != 0;j = nodes[j].father){
					nodes[j].totRound++;
					if (!nodes[j].user)
						nodes[j].winRound++;
				}
				nodes[0].totRound++;
				gameOver = true;
			}
		}

		if (isTie(N,tempTop1)){
			for (int j = sel;j != 0;j = nodes[j].father){
				nodes[j].totRound++;
				if (!nodes[j].user)
					nodes[j].winRound++;
			}
			nodes[0].totRound++;
			gameOver = true;
		}
		// judge finished

		if (gameOver){
			for (int i = 0;i < M;i++)
				for (int j = 0;j < N;j++)
					tempBoard1[i][j] = board[i][j];
			for (int i = 0;i < nn;i++)
				tempTop1[i] = top[i];
			fNode = 0;
		}
		else
			fNode = sel;

		finish = clock();
		totaltime = ((double)finish-start)/CLOCKS_PER_SEC;
	} // end while

	int res = nodes[0].lChild;
	double rate = double(nodes[nodes[0].lChild].winRound)/(nodes[nodes[0].lChild].totRound);
	for (int i = nodes[0].lChild+1;i < nodes[0].rChild;i++){
		double tempRate = double(nodes[i].winRound)/nodes[i].totRound;
		if (tempRate > rate){
			res = i;
			rate = tempRate;
		}
	}
	y = nodes[res].y;
	x = nodes[res].x;

	delete[] tempTop1;
	delete[] tempTop2;
	for (int i = 0;i < M;i++){
		delete[] tempBoard1[i];
		delete[] tempBoard2[i];
	}
	delete[] tempBoard1;
	delete[] tempBoard2;

// end of my part
	clearArray(M, N, board);
	return new Point(x, y);
}

extern "C" __declspec(dllexport) void clearPoint(Point* p){
	delete p;
	return;
}

void clearArray(int M, int N, int** board){
	for(int i = 0; i < M; i++){
		delete[] board[i];
	}
	delete[] board;
}