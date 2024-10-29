from tensorflow.keras.layers import Dropout, Flatten, Conv2D, Input, BatchNormalization, Dense, MaxPooling2D, Activation, AveragePooling2D, GlobalAveragePooling2D
from tensorflow.keras.regularizers import L2
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers import Adam

def feature_extractor(input_shape, l2_rate = 0.01):
    
    l2 = L2(l2_rate)
    input = Input(shape=input_shape)

    # base model
    conv1 = Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(input)
    conv1 = Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(conv1)
    mpool1 = MaxPooling2D()(conv1)

    conv2 = Conv2D(128, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(mpool1)
    conv2 = Conv2D(128, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(conv2)
    mpool2 = MaxPooling2D()(conv2)

    conv3 = Conv2D(256, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(mpool2)
    conv3 = Conv2D(256, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(conv3)
    mpool3 = MaxPooling2D()(conv3)

    conv4 = Conv2D(512, 3, activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(mpool3)
    conv4 = Conv2D(512, 3, activation='relu', padding='same', kernel_initializer="he_normal", bias_initializer="zeros", kernel_regularizer = l2)(conv4)
    mpool4 = GlobalAveragePooling2D()(conv4)

    
    return Model(inputs = input, outputs= mpool4)

def dual_net(input_shape, n_classes_age, n_classes_gender, l2_rate = 0.01):
    
    l2 = L2(l2_rate)
    input = Input(shape=input_shape)

    # base model
    conv1 = Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(input)
    conv1 = Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv1)
    mpool1 = MaxPooling2D()(conv1)

    conv2 = Conv2D(128, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(mpool1)
    conv2 = Conv2D(128, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv2)
    mpool2 = MaxPooling2D()(conv2)

    conv3 = Conv2D(256, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(mpool2)
    conv3 = Conv2D(256, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv3)
    mpool3 = MaxPooling2D()(conv3)

    conv4 = Conv2D(512, 3, activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(mpool3)
    conv4 = Conv2D(512, 3, activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv4)

    # Net 1 for age
    mpool_5 = GlobalAveragePooling2D()(conv4)
    norm_1 = BatchNormalization()(mpool_5)
    net_1_drop1 = Dropout(0.4)(norm_1)
    dense_1 = Dense(128, activation='relu', kernel_initializer="he_normal", kernel_regularizer = l2)(net_1_drop1) 
    dense_2 = Dense(64, activation='relu', kernel_initializer="he_normal", kernel_regularizer = l2)(dense_1)
    drop_1 = Dropout(0.2)(dense_2)
    out1 = Dense(n_classes_age, activation = 'softmax', kernel_initializer = 'he_normal', name='age')(drop_1)

    # Net 2 for gender
    mpool_6 = GlobalAveragePooling2D()(conv4)
    norm_2 = BatchNormalization()(mpool_6)
    net_2_drop1 = Dropout(0.4)(norm_2)
    out2 = Dense(n_classes_gender, activation = 'sigmoid', kernel_initializer = 'he_normal', name='gender')(net_2_drop1)
    return Model(inputs=input, outputs=[out1, out2])

def get_feature_extractor(weights_path):

    model = feature_extractor((224,224, 1), 0.0001)
    model.compile(optimizer = Adam(lr = 0.0001), loss = 'categorical_crossentropy', metrics = ['accuracy'])
    model.load_weights(weights_path, by_name=True)    

    return model

def single_net(input_shape, classes, l2_rate = 0.01):
    
    l2 = L2(l2_rate)
    input = Input(shape=input_shape)

    # base model
    conv1 = Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(input)
    conv1 = Conv2D(64, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv1)
    mpool1 = MaxPooling2D()(conv1)

    conv2 = Conv2D(128, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(mpool1)
    conv2 = Conv2D(128, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv2)
    mpool2 = MaxPooling2D()(conv2)

    conv3 = Conv2D(256, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(mpool2)
    conv3 = Conv2D(256, kernel_size=(3, 3), activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv3)
    mpool3 = MaxPooling2D()(conv3)

    conv4 = Conv2D(512, 3, activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(mpool3)
    conv4 = Conv2D(512, 3, activation='relu', padding='same', kernel_initializer="he_normal", kernel_regularizer = l2)(conv4)
    mpool4 = MaxPooling2D()(conv4)

    # Inference layer
    mpool_5 = GlobalAveragePooling2D()(conv4)
    norm_1 = BatchNormalization()(mpool_5)
    net_1_drop1 = Dropout(0.4)(norm_1)
    dense_1 = Dense(128, activation='relu', kernel_initializer="he_normal", kernel_regularizer = l2)(net_1_drop1) 
    dense_2 = Dense(64, activation='relu', kernel_initializer='he_normal', kernel_regularizer = l2)(dense_1)
    drop_1 = Dropout(0.2)(dense_2)
    out = Dense(classes, activation = 'softmax')(drop_1)

    return Model(inputs = input, outputs = out)