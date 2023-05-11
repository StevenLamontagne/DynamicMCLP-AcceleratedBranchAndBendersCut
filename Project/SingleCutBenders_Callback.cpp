#include "SingleCutBenders_Callback.h"


void SingleCutBenders_Callback::invoke(const IloCplex::Callback::Context& context)
{
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate): 
	{
		double obj_new = context.getCandidateObjective();
		if (obj_new == obj_val) { duplicate_counter += 1; }
		else { obj_val = obj_new; duplicate_counter = 0; }
		if (duplicate_counter >= 5){ 
			cout << "Looping detected at integer node." << endl;
			context.rejectCandidate();
			break; }

		LazyCutCallback(context);
		double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
		if (LazyCutTimes.size() < 512) { LazyCutTimes.push_back(time); }
		break;
	}
	case (IloCplex::Callback::Context::Id::Relaxation):
	{
		double obj_new = context.getRelaxationObjective();
		if (obj_new == obj_val) { duplicate_counter += 1; }
		else { obj_val = obj_new; duplicate_counter = 0; }
		if (duplicate_counter >= 5) { 
			cout << "Looping detected at fractional node." << endl;
			context.exitCutLoop(); 
			break; 
		}


		UserCutCallback(context);
		double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
		if (UserCutTimes.size() < 512) { UserCutTimes.push_back(time); }
		break;
	}
	}
}


void SingleCutBenders_Callback::LazyCutCallback(const IloCplex::Callback::Context& context)
{
	time_t start;
	time(&start);

	IloEnv env = context.getEnv();

	//Get current solution value, x_tilde
	ArrayXXd x_tilde = ArrayXXd::Constant(data.T, data.M_bar, 0.0);
	for (int t = 0; t < data.T; t++) {
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
	}

	//Calculate I_tilde
	vector<VectorXd> I_tilde;
	for (int t = 0; t < data.T; t++) {
		I_tilde.push_back(VectorXd::Constant(data.P[t], 0.0));
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { I_tilde[t] += data.a[t].col(j_bar); }
		}
	}


	AddCuts(context, I_tilde);

}

void SingleCutBenders_Callback::UserCutCallback(const IloCplex::Callback::Context& context)
{
	time_t start;
	time(&start);

	IloEnv env = context.getEnv();
	
	//Get current solution value, x_tilde
	ArrayXXd x_tilde = ArrayXXd::Constant(data.T, data.M_bar, 0.0);
	for (int t = 0; t < data.T; t++) {
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) { x_tilde(t, j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
	}

	//Calculate I_tilde
	vector<VectorXd> I_tilde;
	for (int t = 0; t < data.T; t++) {
		I_tilde.push_back(VectorXd::Constant(data.P[t], 0.0));
		for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > EPS) { I_tilde[t] += x_tilde(t, j_bar) * data.a[t].col(j_bar); }
		}
	}

	AddCuts(context, I_tilde);

}

void SingleCutBenders_Callback::AddCuts(const IloCplex::Callback::Context& context, const vector<VectorXd>& I_tilde)
{
	IloEnv env = context.getEnv();

	IloExpr lhs(env);
	lhs -= theta;
	IloNum covered = 0;

	switch (cut_type)
	{

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B1 cut
	case SINGLE_CUTS::SingleB1:
	{
		for (int t = 0; t < data.T; t++) {
			covered += (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
			VectorXd uncovered = (I_tilde[t].array() < 1).cast<double>();
			for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
				lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
			}
		}
		break;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	//Single B2 cut
	case SINGLE_CUTS::SingleB2:
	{
		for (int t = 0; t < data.T; t++) {
			covered += (I_tilde[t].array() > 1).matrix().cast<double>().dot(data.weights[t]);
			VectorXd uncovered = (I_tilde[t].array() <= 1).cast<double>();
			for (int j_bar = 0; j_bar < data.M_bar; j_bar++) {
				lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered) + data.Ps[t](j_bar)) * x[t][j_bar];
			}
		}
		break;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	//Will probably throw an error earlier, but just in case
	default:
		cout << "Unrecognised cut type" << endl;
		throw;
		break;
	}



	switch (context.getId()) {
	case (IloCplex::Callback::Context::Id::Candidate):
	{
		if (context.getCandidateValue(lhs) < -covered - EPS) {
			context.rejectCandidate(lhs >= -covered);
		}
		lhs.end();
		break;
	}
	case (IloCplex::Callback::Context::Id::Relaxation):
	{
		if (context.getRelaxationValue(lhs) < -covered - EPS) {
			//context.addUserCut(lhs >= -covered, IloCplex::UseCutForce, IloFalse);
			context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
		}
		lhs.end();
		break;
	}
	}
}
