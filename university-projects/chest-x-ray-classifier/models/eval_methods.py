from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

def confusion_matrix_show(y_val_gender, predicted_gender):
    # Confustion matrix for BASE

    # Assuming y_true and y_pred contain the true and predicted labels respectively
    cm = confusion_matrix(y_val_gender, predicted_gender)

    # Create a heatmap visualization of the confusion matrix
    sns.heatmap(cm, annot=True, cmap='Blues', fmt='g')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title('Confusion Matrix')
    plt.show()

def mapping_tool(df):
    # Define the mapping dictionary
    mapping = {
        0: '15_19_F', 1: '15_19_M', 2: '20_24_F', 3: '20_24_M', 4: '25_29_F', 5: '25_29_M',
        6: '30_34_F', 7: '30_34_M', 8: '35_39_F', 9: '35_39_M', 10: '40_44_F', 11: '40_44_M',
        12: '45_49_F', 13: '45_49_M', 14: '50_54_F', 15: '50_54_M', 16: '55_59_F', 17: '55_59_M',
        18: '60_64_F', 19: '60_64_M', 20: '65_69_F', 21: '65_69_M', 22: '70_74_F', 23: '70_74_M',
        24: '75_79_F', 25: '75_79_M', 26: '80_84_F', 27: '80_84_M', 28: '85_89_F', 29: '85_89_M'
    }

    # Map the values using the mapping dictionary
    df['True'] = df['True'].map(mapping)
    df['Predicted'] = df['Predicted'].map(mapping)
    
    return df

def graph_heat(df):
  # Calculate the frequency of each point
  point_frequency = df.groupby(['True', 'Predicted']).size().reset_index(name='Frequency')

  # Create a pivot table to reshape the data for the heatmap
  heatmap_data = point_frequency.pivot('Predicted', 'True', 'Frequency')

  # Set up the heatmap figure
  fig, ax = plt.subplots(figsize=(20, 10))

  # Invert the 'Predicted' axis in the heatmap data
  heatmap_data = heatmap_data.iloc[::-1]

  # Create the heatmap
  sns.heatmap(heatmap_data, cmap='YlGnBu', annot=True, cbar=True, ax=ax, fmt='.0f')

  # Add labels and title
  plt.xlabel('True Age Cohort and Gender')
  plt.ylabel('Predicted Age Cohort and Gender')
  plt.title('Heatmap based on Frequency')

  # Show the plot
  plt.show()