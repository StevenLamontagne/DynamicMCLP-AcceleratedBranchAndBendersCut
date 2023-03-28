#include "Callback_LocalBranching.h"



void Callback_LocalBranching::invoke(const IloCplex::Callback::Context& context)
{
	try {

		double BestBound = context.getDoubleInfo(IloCplex::Callback::Context::Info::BestBound);
		if (BestBound < UpperBound) {
			cout << "Dual bound: " << BestBound << endl;
			UpperBound = BestBound;
			if (UpperBound < thresholdObjective) {
				cout << "Optimality condition reached." << endl;
				status = CallbackStatus::Optimal;
				context.abort();
			}
		}

		switch (context.getId()) {
		case (IloCplex::Callback::Context::Id::Candidate):
		{
			double obj = context.getCandidateObjective();

			if (obj < LowerBound) {
				context.rejectCandidate();
				cout << "Candidate solution rejected." << endl;
				return;
			}


			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t, j_bar) = context.getCandidatePoint(x[t][j_bar]); }
				theta_tilde(t) = context.getCandidatePoint(theta[t]);
			}
			x_tilde = x_tilde.round();

			//if (obj < LowerBound) {
			//	IloEnv env = context.getEnv();
			//	IloExpr taboo(env);
			//	for (int t = 0; t < T; t++) {
			//		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			//			if (x_tilde(t, j_bar) > 1 - EPS) { taboo += 1 - x[t][j_bar]; }
			//			else { taboo += x[t][j_bar]; }
			//		}
			//	}
			//	context.rejectCandidate(taboo >= 1);
			//	cout << "Candidate solution rejected." << endl;
			//	return;
			//}
			cout << "Candidate objective bound: " << obj << endl;
			//Get current solution value, x_tilde


			status = CallbackStatus::Integer;
			context.abort();
			break;
		}
		case (IloCplex::Callback::Context::Id::Relaxation):
		{
			
			double obj = context.getRelaxationObjective();
			if (obj < LowerBound) { 
				cout << "Fractional solution rejected." << endl;
				context.pruneCurrentNode();
				return;
			}
			//cout << "Fractional objective bound: " << obj_new << endl;
			//if ((SkipFractional) || (context.getIntInfo(IloCplex::Callback::Context::Info::NodeUID) != 0)) { break; }

			//chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
			//IloEnv env = context.getEnv();

			//Get current solution value, x_tilde
			for (int t = 0; t < T; t++) {
				for (int j_bar = 0; j_bar < M_bar; j_bar++) { x_tilde(t,j_bar) = context.getRelaxationPoint(x[t][j_bar]); }
			}
			
			if (WithinRestrictedDistanceFast(x_tilde, x_tilde.floor())) {
				x_tilde = x_tilde.floor();
				for (int t = 0; t < T; t++) {
					theta_tilde(t) = context.getRelaxationPoint(theta[t]);
				}

				status = CallbackStatus::FractionalFloor;
				context.abort();
			}
			else if (WithinRestrictedDistance(x_tilde, x_tilde.round())) {
				x_tilde = x_tilde.round();
				for (int t = 0; t < T; t++) {
					theta_tilde(t) = context.getRelaxationPoint(theta[t]);
				}
				status = CallbackStatus::FractionalRound;
				context.abort();
			}
			//FractionalCutCallback(context, x_tilde);
			//stats["nUserCuts"] += 1;
			//double time = chrono::duration_cast<chrono::duration<double>>(chrono::steady_clock::now() - t1).count();
			//if (UserCutTimes.size() < 512) { UserCutTimes.push_back(time); }
			//cout << "User cut time: " <<  time << endl;

			break;
		}
		case (IloCplex::Callback::Context::Id::Branching):
		{			
			double obj = context.getRelaxationObjective();
			if (obj < LowerBound) {
				cout << "Branch pruned." << endl;
				context.pruneCurrentNode();
				return;
		}
		}
		default:
			return;
		}
	}
	catch (IloException & e) {
		cout << "Exception during callback: " << e << endl;
		throw e;
	}
	//cout << "Returning to master problem." << endl;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Callback_LocalBranching::FractionalCutCallback(const IloCplex::Callback::Context& context, ArrayXXd& x_tilde)
{
	//IloEnv env = context.getEnv();
	//if (WithinRestrictedDistanceFast(x_tilde, x_tilde.floor())) {
	//	x_tilde = x_tilde.floor();
	//	for (int t = 0; t < T; t++) {
	//		theta_tilde(t) = context.getRelaxationPoint(theta[t]);
	//	}
	//	
	//	status = CallbackStatus::FractionalFloor;
	//	context.abort();
	//}
	//else if (WithinRestrictedDistance(x_tilde, x_tilde.round())) {
	//	x_tilde = x_tilde.round();
	//	for (int t = 0; t < T; t++) {
	//		theta_tilde(t) = context.getRelaxationPoint(theta[t]);
	//	}
	//	status = CallbackStatus::FractionalRound;
	//	context.abort();
	//}
	////else{ AddBendersCuts(context, x_tilde); }
	AddBendersCuts(context, x_tilde); 
}

void Callback_LocalBranching::AddBendersCuts(const IloCplex::Callback::Context& context, const ArrayXXd& x_tilde)
{
	IloEnv env = context.getEnv();
	vector<VectorXd> I_tilde = CalculateItilde(x_tilde);

	////Multi-cut (by year), B1
//for (int t = 0; t < T; t++) {
//	IloExpr lhs(env);
//	lhs -= theta[t];
//	double covered = (I_tilde[t].array() >= 1).matrix().cast<double>().dot(data.weights[t]);
//	VectorXd uncovered = (I_tilde[t].array() < 1).matrix().cast<double>();
//	for (int j_bar = 0; j_bar < M_bar; j_bar++) {
//		lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
//	}

//	switch (context.getId()) {
//	case (IloCplex::Callback::Context::Id::Candidate):
//		if (context.getCandidateValue(lhs) < -covered -EPS) {
//			//nCutsAdded += 1;
//			context.rejectCandidate(lhs >= -covered);
//		}
//		lhs.end();
//		break;
//	case (IloCplex::Callback::Context::Id::Relaxation):
//		if (context.getRelaxationValue(lhs) < -covered -EPS) {
//			//context.addUserCut(lhs >= 0, IloCplex::UseCutForce, IloFalse).end();
//			context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
//		}
//		lhs.end();
//		break;
//	}
//}



//Multi-cut (by year), Multi1PO1
	for (int t = 0; t < T; t++) {
		IloExpr lhs(env);
		lhs -= theta[t];
		double covered = ((I_tilde[t].array() > 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() >= 1))).matrix().cast<double>().dot(data.weights[t]);

		VectorXd uncovered = ((I_tilde[t].array() < 1) || ((I_tilde[t].array() == 1) && (core_coverage[t].array() < 1))).cast<double>();
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			lhs += (data.CutCoeffs[t].col(j_bar).dot(uncovered)) * x[t][j_bar];
		}

		if (context.getRelaxationValue(lhs) < -covered - EPS) {
			//context.addUserCut(lhs >= -covered, IloCplex::UseCutForce, IloFalse);
			context.addUserCut(lhs >= -covered, IloCplex::UseCutPurge, IloFalse);
		}
		lhs.end();
	}

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//Calculate the l-infinity norm across all years for each j_bar, then return the sum 
bool Callback_LocalBranching::WithinRestrictedDistance(const ArrayXXd& x1, const ArrayXXd& x2)
{
	for (int t = 0; t < T; t++) {
		if ((x1 - x2).abs().sum() > 2.0) { return false; }
	}
	return true;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//As the CalculateDistance function, but when we are guaranteed that x_upper >= x 
bool Callback_LocalBranching::WithinRestrictedDistanceFast(const ArrayXXd& x_upper, const ArrayXXd& x)
{
	for (int t = 0; t < T; t++) {
		if ((x_upper - x).sum() > 2.0) { return false; }
	}
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



double Callback_LocalBranching::CalculateObjective(const ArrayXXd& x_tilde)
{
	double covered = 0.0;
	for (int t = 0; t < T; t++) {
		covered += data.Precovered[t];
		VectorXd I_tilde = VectorXd::Zero(P[t]);
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { 
				I_tilde += data.a[t].col(j_bar); 
				covered += data.Ps[t](j_bar);
			}
		}
		covered += (I_tilde.array() >= 1).matrix().cast<double>().dot(data.weights[t]);
	}

	return covered;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<VectorXd> Callback_LocalBranching::CalculateItilde(const ArrayXXd& x_tilde)
{
	vector<VectorXd> I_tilde;
	for (int t = 0; t < T; t++) {
		I_tilde.push_back(VectorXd::Constant(P[t], 0.0));
		for (int j_bar = 0; j_bar < M_bar; j_bar++) {
			if (x_tilde(t, j_bar) > 1 - EPS) { I_tilde[t] += data.a[t].col(j_bar); }
		}
	}
	return I_tilde;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

