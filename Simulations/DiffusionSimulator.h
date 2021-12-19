#ifndef DIFFUSIONSIMULATOR_h
#define DIFFUSIONSIMULATOR_h

#include "Simulator.h"
#include "vectorbase.h"

//impement your own grid class for saving grid data
class Grid {
public:
	// Construtors
	Grid(size_t n, size_t m);

	void reset() {
		for (size_t i = 0; i < N; i++)
			for (size_t j = 0; j < M; j++)
				set(i,j,rand()%2*200-100);
	}
	
	void set(int x, int y, double value) {
		if (x == 0 || y == 0 || x == N - 1 || y == M - 1) return;
		matrix[x][y] = value;
	}


	double operator()(size_t i, size_t j) {
		return matrix[i][j];
	}

	double operator()(size_t i) {
		const int row = i / M;
		const int col = i % M;
		return matrix[row][col];
	}

	int getN() {
		return N;
	}
	int getM() {
		return M;
	}

private:
	// Attributes
	std::vector<std::vector<double>> matrix;
	int M, N;
};



class DiffusionSimulator:public Simulator{
public:
	// Construtors
	DiffusionSimulator();

	// Functions
	const char * getTestCasesStr();
	void initUI(DrawingUtilitiesClass * DUC);
	void reset();
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext);
	void notifyCaseChanged(int testCase);
	void simulateTimestep(float timeStep);
	void externalForcesCalculations(float timeElapsed) {};
	void onClick(int x, int y);
	void onMouse(int x, int y);
	// Specific Functions
	void drawObjects();
	void diffuseTemperatureExplicit(float factor);
	void diffuseTemperatureImplicit(float timestep);

private:
	// Attributes
	Vec3  m_vfMovableObjectPos;
	Vec3  m_vfMovableObjectFinalPos;
	Vec3  m_vfRotate;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
	//Grid *T; //save results of every time step
	Grid T;
	int M, N;

	double alpha;
};

#endif