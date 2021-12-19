#include "DiffusionSimulator.h"
#include "pcgsolver.h"
using namespace std;

Grid::Grid(size_t n, size_t m) : N(n), M(m) {
	matrix.resize(n);
	for (size_t i = 0; i < n; i++) matrix[i].resize(m);
}

DiffusionSimulator::DiffusionSimulator() : M(11), N(11), T(11,11), alpha(0.1)
{
	m_iTestCase = 0;
	m_vfMovableObjectPos = Vec3();
	m_vfMovableObjectFinalPos = Vec3();
	m_vfRotate = Vec3();
}

const char * DiffusionSimulator::getTestCasesStr(){
	return "Explicit_solver, Implicit_solver";
}

void DiffusionSimulator::reset(){
	m_mouse.x = m_mouse.y = 0;
	m_trackmouse.x = m_trackmouse.y = 0;
	m_oldtrackmouse.x = m_oldtrackmouse.y = 0;

	T.reset();
	T.set(T.getN() / 2, T.getM() / 2, 100);
}

void DiffusionSimulator::initUI(DrawingUtilitiesClass * DUC)
{
	this->DUC = DUC;
	TwAddVarRW(DUC->g_pTweakBar, "alpha", TW_TYPE_DOUBLE, &alpha, "step=0.05 min=0.05");
	TwAddVarRW(DUC->g_pTweakBar, "N", TW_TYPE_INT32, &N, "step=1 min=5");
	TwAddVarRW(DUC->g_pTweakBar, "M", TW_TYPE_INT32, &M, "step=1 min=5");

	/*switch (m_iTestCase)
	{
	case 0:break;
	case 1:
		break;
	case 2:break;
	default:break;
	}*/
}

void DiffusionSimulator::notifyCaseChanged(int testCase)
{
	m_iTestCase = testCase;
	m_vfMovableObjectPos = Vec3(0, 0, 0);
	m_vfRotate = Vec3(0, 0, 0);
	//
	//to be implemented
	//
	reset();
	switch (m_iTestCase)
	{
	case 0:
		cout << "Explicit solver!\n";
		break;
	case 1:
		cout << "Implicit solver!\n";
		break;
	default:
		cout << "Empty Test!\n";
		break;
	}
}

void DiffusionSimulator::diffuseTemperatureExplicit() {//add your own parameters

	// to be implemented
	//make sure that the temperature in boundary cells stays zero
	
}

void setupB(std::vector<Real>& b, Grid& T) {//add your own parameters
	// to be implemented
	//set vector B[sizeX*sizeY]
	const int S = T.getN() * T.getM();
	for (int i = 0; i < S; i++) {
		b.at(i) = T(i);
	}
}

void fillT(std::vector<Real> x, Grid& T) {//add your own parameters
	// to be implemented
	//fill T with solved vector x
	//make sure that the temperature in boundary cells stays zero
	for (size_t i = 0; i < T.getN(); i++)
		for (size_t j = 0; j < T.getM(); j++)
			T.set(i, j, x[i * T.getM() + j]);	
}

void setupA(SparseMatrix<Real>& A, double factor, int N, int M) {//add your own parameters
	// to be implemented
	//setup Matrix A[sizeX*sizeY*sizeZ, sizeX*sizeY*sizeZ]
	// set with:  A.set_element( index1, index2 , value );
	// if needed, read with: A(index1, index2);
	// avoid zero rows in A -> set the diagonal value for boundary cells to 1.0

	const int S = N * M;
	for (int i = 0; i < S; i++) {
		A.set_element(i, i, 1); // set diagonal
	}

	for (size_t i = 1; i < N-1; i++) {
		for (size_t j = 1; j < M-1; j++) {
			const size_t ind = i * M + j;
			A.set_element(ind, ind - M, -factor);
			A.set_element(ind, ind - 1, -factor);
			A.set_element(ind, ind, 4 * factor + 1);
			A.set_element(ind, ind + 1, -factor);
			A.set_element(ind, ind + M, -factor);
		}
	}
}


void DiffusionSimulator::diffuseTemperatureImplicit(float timestep) {//add your own parameters
	// solve A T = b
	// to be implemented
	const int S = T.getM() * T.getN(); //N = sizeX*sizeY*sizeZ
	SparseMatrix<Real> A(S);
	std::vector<Real> b(S);
	setupA(A, timestep * alpha, T.getN(), T.getM());
	setupB(b, T);

	// perform solve
	Real pcg_target_residual = 1e-05;
	Real pcg_max_iterations = 1000;
	Real ret_pcg_residual = 1e10;
	int  ret_pcg_iterations = -1;

	SparsePCGSolver<Real> solver;
	solver.set_solver_parameters(pcg_target_residual, pcg_max_iterations, 0.97, 0.25);

	std::vector<Real> x(S);
	for (int j = 0; j < S; ++j) { x[j] = 0.; }

	// preconditioners: 0 off, 1 diagonal, 2 incomplete cholesky
	solver.solve(A, b, x, ret_pcg_residual, ret_pcg_iterations, 0);

	// x contains the new temperature values
	fillT(x, T);//copy x to T
}



void DiffusionSimulator::simulateTimestep(float timeStep)
{
	// to be implemented
	// update current setup for each frame

	if (N != T.getN() || M != T.getM()) {
		T = Grid(N, M);
		reset();
	}

	switch (m_iTestCase)
	{
	case 0:
		//delete T;
		//T = diffuseTemperatureExplicit();
		break;
	case 1:
		diffuseTemperatureImplicit(timeStep);
		break;
	}
}

void DiffusionSimulator::drawObjects()
{
	// to be implemented
	//visualization
	const double cX = T.getN()/2;
	const double cY = T.getM()/2;
	for (size_t i = 0; i < T.getN(); i++) {
		for (size_t j = 0; j < T.getM(); j++) {
			const double val = T(i, j);
			const double red = (val + 100) / 200;
			const double blue = 1 - red;
			if (i == 0 || j == 0 || i == T.getN() - 1 || j == T.getM() - 1)
				DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, Vec3(0, 0, 0));
			else
				DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, Vec3(red, 0, blue));
			DUC->drawSphere(Vec3((-cX+i)/ T.getN(), (-cY+j)/ T.getM(), 0), Vec3(1, 1, 1) * 0.05f);
		}
	}
}


void DiffusionSimulator::drawFrame(ID3D11DeviceContext* pd3dImmediateContext)
{
	drawObjects();
}

void DiffusionSimulator::onClick(int x, int y)
{
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}

void DiffusionSimulator::onMouse(int x, int y)
{
	m_oldtrackmouse.x = x;
	m_oldtrackmouse.y = y;
	m_trackmouse.x = x;
	m_trackmouse.y = y;
}
