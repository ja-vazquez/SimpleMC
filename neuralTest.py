from simplemc.DriverMC import DriverMC
#import tensorflow as tf

# inputs = tf.keras.Input(shape=(3,))
# x = tf.keras.layers.Dense(300, activation=tf.nn.relu)(inputs)
# outputs = tf.keras.layers.Dense(1, activation=tf.nn.relu)(x)
# model = tf.keras.Model(inputs=inputs, outputs=outputs)

analyzer = DriverMC(analyzername="nested", model="LCDM", datasets="HD", chainsdir="chains")

analyzer.nestedRunner(neuralNetwork=True, nlivepoints=50)
analyzer.postprocess()

# analyzer = DriverMC(iniFile='baseConfig.ini')
# analyzer.executer()
# analyzer.postprocess()