import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, MinMaxScaler

def train_model(training_set_path):
    """
    Train a RandomForestClassifier using the provided training set and return the trained model.
    
    Args:
        training_set_path (str): Path to the training set CSV file.
        
    Returns:
        RandomForestClassifier: Trained RandomForestClassifier model.
    """
    # Read the training set CSV file into a pandas dataframe
    acr_df = pd.read_csv(training_set_path)

    # Define input features and output class
    X = acr_df.drop(['Name', 'DP'], axis=1).values
    y = acr_df['DP'].values

    # Define the machine learning model pipeline
    pipelineDP = Pipeline(steps=[
        ('PFeatures', PolynomialFeatures(2)),
        ('scaler', MinMaxScaler()),
        ('model', RandomForestClassifier(n_estimators=1500, criterion="entropy", max_features='auto',
                                         max_depth=30, bootstrap=True, min_samples_leaf=3,
                                         oob_score=False, min_samples_split=2))
    ])

    # Fit the machine learning model
    model_rf = pipelineDP.fit(X, y)

    return model_rf

def predict_depolymerases(protein_df, model):
    """
    Predict depolymerase probabilities for each protein in the dataframe.
    
    Args:
        protein_df (pd.DataFrame): DataFrame containing protein data.
        model (RandomForestClassifier): Trained RandomForestClassifier model.
        
    Returns:
        pd.DataFrame: DataFrame containing depolymerase probabilities for each protein.
    """
    # Predict probabilities of DePol for each protein
    probabilities = model.predict_proba(protein_df.drop('name', axis=1))

    # Create a new dataframe with the results
    results_df = protein_df[['name']].copy()
    results_df['Probability_DePol'] = probabilities[:, 1]

    return results_df