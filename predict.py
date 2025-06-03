import pandas as pd
import lightgbm as lgb

print("📥 正在加载模型...")
model = lgb.Booster(model_file="SaGP.model")

print("📄 正在读取 ACC_PSSM.csv 文件...")
df_new = pd.read_csv("ACC_PSSM.csv")
df_new.columns = df_new.columns.str.replace(r"[^\w]", "_", regex=True)

acc_pssm_features = [col for col in df_new.columns if "ACC_PSSM" in col]
X_new = df_new[acc_pssm_features]
print(f"🔢 样本数: {X_new.shape[0]}, 特征数: {X_new.shape[1]}")

print("🤖 正在进行预测...")
y_pred_prob = model.predict(X_new)
y_pred = (y_pred_prob >= 0.5).astype(int)

output = pd.DataFrame({
    "Sample_ID": df_new.index,
    "Pred_Probability": y_pred_prob,
    "Pred_Label": y_pred
})

output.to_csv("prediction_results.csv", index=False)
print("✅ 预测完成，结果已保存至 prediction_results.csv")
