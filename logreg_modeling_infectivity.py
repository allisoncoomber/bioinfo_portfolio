# infectivity - proportion of inoculated plants which became infected
# this is a better way to describe what I have been calling extinctions

# factors
# H = host genotype
# L = pathogen lineage
# R = replication
# u = mean
# e = error

## was there an infection = u + H + L + H*L + e

import pandas as pd
import sklearn as sklearn
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

extinctions = pd.read_csv(r"path")
hosts = pd.get_dummies(extinctions.Host, prefix='Host')
pathogens = pd.get_dummies(extinctions.Pathogen, prefix='Pathogen')
treatments = pd.get_dummies(extinctions.Treatment, prefix='Treatment')

#expand categorial columns into dummy variables
extinctions=extinctions.join(hosts)
extinctions=extinctions.join(pathogens)
extinctions=extinctions.join(treatments)
extinctions.drop(['Host','Pathogen','Treatment'], axis=1, inplace=True)


X=extinctions.drop(['Extinctions'], axis=1)
y=extinctions["Extinctions"].copy()

X_train,X_test,y_train,y_test=train_test_split(X,y,test_size=0.25,random_state=0)

#instantiate model
logreg=LogisticRegression()
#fit model
logreg.fit(X_train,y_train)

print(list(zip(logreg.coef_[0], X_train.columns)))

#predict outcomes from test set
y_pred=pd.Series(logreg.predict(X_test))
# check out our predictions
y_test = y_test.reset_index(drop=True)
z = pd.concat([y_test, y_pred], axis=1)
z.columns = ['True', 'Prediction']

#calculate feature importance scores for different independent variables
# get importance
importance = logreg.coef_[0]
# summarize feature importance
#for i,v in enumerate(importance):
#	print('Feature: %0d, Score: %.5f' % (i,v))
# plot feature importance
fig, ax = plt.subplots()
ax.set_xticks(np.arange(36), labels=X_train.columns)
plt.xticks(rotation = 90)
fig.subplots_adjust(bottom=0.65)
plt.bar([x for x in range(len(importance))], importance)
plt.show()


#calculate some stats on the model
print("Accuracy:", metrics.accuracy_score(y_test, y_pred))
print("Precision:", metrics.precision_score(y_test, y_pred))
print("Recall:", metrics.recall_score(y_test, y_pred))


#make a confusion matrix
cnf_matrix = metrics.confusion_matrix(y_test, y_pred)

labels = [0, 1]
fig, ax = plt.subplots()
tick_marks = np.arange(len(labels))
plt.xticks(tick_marks, labels)
plt.yticks(tick_marks, labels)
# create heatmap
sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu", fmt='g')
ax.xaxis.set_label_position("top")
plt.title('Confusion matrix', y=1.1)
plt.ylabel('True')
plt.xlabel('Predicted')
plt.show()