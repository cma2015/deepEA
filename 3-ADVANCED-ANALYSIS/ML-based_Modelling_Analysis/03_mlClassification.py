import matplotlib, math,pickle, argparse
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc, precision_recall_curve,average_precision_score
from scipy import interp
from xgboost import XGBClassifier

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-positive", dest = "posMat", type = str, default = None,
                        help = "The feature matrix of positive samples.")
    parser.add_argument("-negative", dest ="negMat", type = str, default = None,
                        help = "The feature matrix of negative samples.")
    parser.add_argument("-m", dest = "method", help = "Machine learning algorithms.")
    parser.add_argument("-kernel", dest = "kernel", help = "Kernel function used for SVM.")
    parser.add_argument("-gamma", dest = "gamma", type=str, default = "None",
                        help = "Kernel coefficient for rbf, poly and sigmoid.")
    parser.add_argument("-coef0", dest = "coef0", type=float, default = 0.0,
                        help = "Independent term in kernel function. It is only significant in poly and sigmoid.")
    parser.add_argument("-degree", dest = "degree", type=int, default=3,
                        help = "Degree of the polynomial kernel function poly. Ignored by all other kernels.")
    parser.add_argument("-lr", dest = "lr", type=float, default=0.01,
                        help = "Learning rate required for XGBoost classifier.")
    parser.add_argument("-solver", dest = "solver", type = str, default = "liblinear",
                        help = "Algorithm to use in optimization problem.")
    parser.add_argument("-out", dest = "outDir", type = str, default = None,
                        help = "The output directory.")
    parser.add_argument("-t", dest = "cpus", type = int, default = 1,
                        help="The number of threads.")
    parser.add_argument("-k", dest="k", type=int, default = 0,
                        help = "The k-fold cross validation. Default: 0 (not perform cross-validation).")
    parser.add_argument("-threshold", dest="threshold", type=float, default = None,
                        help = "The threshold used for classifying CMRs.")
    parser.add_argument("-estimators", dest="n_estimators", type=int, default = 500,
                        help = "The number of trees used for build random forest classifier.")
    args = parser.parse_args()
    return args


def cross_validation(X, y, k, clf):
    cv = StratifiedKFold(n_splits = k)
    res = {}
    i=1
    for train, test in cv.split(X, y):
        clf.fit(X[train], y[train])
        yscore = clf.predict_proba(X[test])
        tmpID = "fold_" + str(i)
        curDic = {}
        curDic["yscore"] = yscore
        curDic["ytest"] = y[test]
        res[tmpID] = curDic
        i = i + 1
    return res

def evalModel(posScore, negScore, threshold = 0.5, beta = 1):
    TP = float(sum(posScore > threshold))
    TN = float(sum(negScore <= threshold))
    FN = float(len(posScore)-TP)
    FP = float(len(negScore)-TN)
    res = {}
    res['Sn'] = TP/(TP + FN)
    res['Sp'] = TN/(TN + FP)
    res['Pr'] = TP/(TP + FP)
    res['Acc'] = (TP+TN)/(TP+TN+FP+FN)
    res['Fscore'] = ((1+beta*beta)*res['Pr']*res['Sn'])/(beta*beta*res['Pr']+res['Sn'])
    res['MCC']=(TP*TN-FP*FN)/math.sqrt(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    return res

def visualizeCV(cvRes, beta = 1):
    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3)
    ax4 = fig.add_subplot(2, 2, 4)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    i = 1
    for fold in cvRes.keys():
        curFold = cvRes[fold]
        yscore = curFold["yscore"][:, 1]
        ytest = curFold["ytest"]
        fpr, tpr, thresholds = roc_curve(ytest, yscore)
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        ax1.plot(fpr, tpr, lw=1, alpha=0.3,
                 label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
        i += 1
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax1.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax1.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')
    ax1.set_xlim([-0.05, 1.05])
    ax1.set_ylim([-0.05, 1.05])
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title('Receiver operating characteristic example')
    ax1.legend(loc="lower right")

    yscores = []
    ytests = []
    i = 1
    for fold in cvRes.keys():
        curFold = cvRes[fold]
        yscore = curFold["yscore"][:, 1]
        ytest = curFold["ytest"]
        precision, recall, _ = precision_recall_curve(ytest, yscore)
        yscores.append(yscore)
        ytests.append(ytest)
        curAUC = average_precision_score(ytest, yscore)
        ax2.plot(recall, precision, lw=1, alpha=0.3, label='PR fold %d (prAUC = %0.2f)' % (i, curAUC))
        i += 1

    yscores = np.concatenate(yscores)
    ytests = np.concatenate(ytests)
    precision, recall, _ = precision_recall_curve(ytests, yscores)
    ax2.plot(recall, precision, color='b', lw=2, alpha=.8)
    ax2.set_xlim([-0.05, 1.05])
    ax2.set_ylim([-0.05, 1.05])
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_title('precision-recall curve')
    ax2.legend(loc="lower right")

    thresVec = [float(i + 1)/float(1000) for i in range(900)]
    thresVec = thresVec[199:(len(thresVec) - 100)]
    posScores = []
    negScores = []
    for fold in cvRes.keys():
        curFold = cvRes[fold]
        yscore = curFold["yscore"][:, 1]
        ytest = curFold["ytest"]
        pos_score = list(yscore[np.where(ytest == 1)])
        neg_score = list(yscore[np.where(ytest == 0)])
        posScores.extend(pos_score)
        negScores.extend(neg_score)
    posScores = np.array(posScores)
    negScores = np.array(negScores)
    FscoreVec = []
    for j in thresVec:
        curFscore = evalModel(posScore=posScores, negScore=negScores, threshold=j, beta=beta)["Fscore"]
        FscoreVec.append(curFscore)
    ax3.plot(thresVec, FscoreVec, lw=1, alpha=0.3)
    maxF = max(FscoreVec)
    threshold = thresVec[FscoreVec.index(maxF)]
    ax3.set_title("Threshold:" + str(threshold))
    ax3.set_xlabel('Threshold')
    ax3.set_ylabel('F1-score')

    thresVec = [float(i + 1)/float(1000) for i in range(1000)]
    thresVec = thresVec[199:(len(thresVec) - 100)]
    posScores = []
    negScores = []
    for fold in cvRes.keys():
        curFold = cvRes[fold]
        yscore = curFold["yscore"][:, 1]
        ytest = curFold["ytest"]
        pos_score = list(yscore[np.where(ytest == 1)])
        neg_score = list(yscore[np.where(ytest == 0)])
        posScores.extend(pos_score)
        negScores.extend(neg_score)
    posScores = np.array(posScores)
    negScores = np.array(negScores)
    Sn = evalModel(posScore=posScores, negScore=negScores, threshold=threshold)["Sn"]
    Sp = evalModel(posScore=posScores, negScore=negScores, threshold=threshold)["Sp"]
    Pr = evalModel(posScore=posScores, negScore=negScores, threshold=threshold)["Pr"]
    Acc = evalModel(posScore=posScores, negScore=negScores, threshold=threshold)["Acc"]
    MCC = evalModel(posScore=posScores, negScore=negScores, threshold=threshold)["MCC"]
    Fscore = evalModel(posScore=posScores, negScore=negScores, threshold=threshold)["Fscore"]
    x = ['Sn', 'Sp', 'Pr', 'Acc', 'MCC', 'Fscore']
    y = [Sn, Sp, Pr, Acc, MCC, Fscore]
    width = 0.75  # the width of the bars
    ind = np.arange(len(y))  # the x locations for the groups
    ax4.barh(ind, y, width, color=['#FF3333', '#FFCC33', '#99CC33', '#66CC33', '#66CCCC', '#336699'])
    for i, v in enumerate(y):
        ax4.text(v + 3, i + .25, str(v), color='black', fontweight='bold')
    ax4.set_yticks(ind + width / 2)
    ax4.set_yticklabels(x, minor=False)
    ax4.set_title('Evaluation in test samples.')
    ax4.set_xlabel("")
    ax4.set_ylabel('Measure')
    return fig

if __name__ == '__main__':
    args = parse_args()
    outDir = args.outDir
    posFeatureMat = np.loadtxt(args.posMat, dtype = "str")
    negFeatureMat = np.loadtxt(args.negMat, dtype = "str")
    featureName = posFeatureMat[0]
    posFeatureMat = np.delete(posFeatureMat, 0, axis=0)
    negFeatureMat = np.delete(negFeatureMat, 0, axis=0)
    featureMat = np.concatenate((posFeatureMat, negFeatureMat), axis = 0)
    featureMat = featureMat.astype(np.float64)
    label = np.array([1]*posFeatureMat.shape[0] + [0]*negFeatureMat.shape[0])

    # build a classifier
    method = args.method
    if method == "randomForest":
        n_estimators = args.n_estimators
        cpus = args.cpus
        clf = RandomForestClassifier(n_estimators = n_estimators,
                                    criterion = "gini",
                                    n_jobs = cpus)
    elif method == "SVM":
        kernel = args.kernel
        if args.gamma == "None":
            gamma = "auto"
        else:
            gamma = args.gamma
        if kernel == "poly":
            clf = SVC(kernel = "poly", degree = args.degree,
                    coef0 = args.coef0, gamma = gamma, probability = True)
        elif kernel == "rbf" or kernel == "sigmoid":
            clf = SVC(kernel = kernel, gamma = gamma, probability = True)
        else:
            clf = SVC(kernel = kernel, probability = True)
    elif method == "XGBoost":
        clf = XGBClassifier(random_state = 1, learning_rate = args.lr, n_jobs = args.cpus, n_estimators = args.n_estimators)
    elif method == "LogisticRegression":
        clf = LogisticRegression(random_state = 0, solver = args.solver, n_jobs = args.cpus)
    else:
        clf = DecisionTreeClassifier()

    if args.k != 0:
        cvRes = cross_validation(X = featureMat, y = label, k = args.k, clf = clf)
        fig = visualizeCV(cvRes = cvRes, beta = 1)
        fig.savefig(outDir + "CV_evaluation.pdf")
        clf = clf.fit(featureMat, label)
    else:
        clf = clf.fit(featureMat, label)

    with open(outDir + "model.pkl", "wb") as f:
        pickle.dump(clf, f, pickle.HIGHEST_PROTOCOL)   

