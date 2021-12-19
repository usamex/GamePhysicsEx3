#ifndef DIFFUSIONSIMULATOR_h
#define DIFFUSIONSIMULATOR_h

#include "Simulator.h"
#include "vectorbase.h"

//impement your own grid class for saving grid data
class Grid {
public:
	// Construtors
	Grid();
	void changeMatrixElement(int x, int y, double value) {
		matrix.at(x).at(y) = value;
	}
	void setMatrixValue(double value) {
		for (size_t i = 0; i < matrix.size(); i++)
		{
			std::fill(matrix.at(i).begin(), matrix.at(i).end(), value);
		}
	}
	void setBoundaryValues() {
		std::fill(matrix.at(0).begin(), matrix.at(0).end(), 0);
		std::fill(matrix.at(matrix.size() - 1).begin(), matrix.at(matrix.size() - 1).end(), 0);

		for (size_t i = 1; i < matrix.size() - 1; i++)
		{
			matrix.at(i).at(0) = 0;
			matrix.at(i).at(matrix.at(i).size() - 1) = 0;
		}
	}

	int getN() {
		return N;
	}
	int getM() {
		return M;
	}
	void setN(int n) {
		N = n;
	}
	void setM(int m) {
		M = m;
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
	Grid* diffuseTemperatureExplicit();
	void diffuseTemperatureImplicit();

private:
	// Attributes
	Vec3  m_vfMovableObjectPos;
	Vec3  m_vfMovableObjectFinalPos;
	Vec3  m_vfRotate;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
	//Grid *T; //save results of every time step
	Grid *T;
	int M, N;
};

#endif