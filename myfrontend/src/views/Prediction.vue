<template>
  <div class="prediction-container">
    <div class="input-section">
      <el-form :model="predictionForm" label-width="120px">
        <el-form-item label="溶剂SMILES">
          <el-input
            v-model="predictionForm.solventSmiles"
            placeholder="请输入溶剂的SMILES分子式"
            :disabled="!isLoggedIn"
          />
        </el-form-item>
        <el-form-item label="溶质SMILES">
          <el-input
            v-model="predictionForm.soluteSmiles"
            placeholder="请输入溶质的SMILES分子式"
            :disabled="!isLoggedIn"
          />
        </el-form-item>
        <el-form-item>
          <el-button type="primary" @click="handlePredict" :disabled="!isLoggedIn">
            开始预测
          </el-button>
        </el-form-item>
      </el-form>
    </div>

    <div class="result-section" v-if="predictionResult">
      <el-card class="result-card">
        <template #header>
          <div class="card-header">
            <span>预测结果</span>
          </div>
        </template>
        <div class="result-content">
          <div class="result-item">
            <span class="label">相容性评分：</span>
            <span class="value">{{ predictionResult.score }}</span>
          </div>
          <div class="result-item">
            <span class="label">预测说明：</span>
            <p class="value">{{ predictionResult.explanation }}</p>
          </div>
          <div class="result-item">
            <span class="label">分子结构：</span>
            <div class="molecule-images">
              <div class="molecule-image-container">
                <img :src="predictionResult.solventImage" alt="溶剂结构" class="molecule-image" />
                <div class="image-label">溶剂结构</div>
              </div>
              <div class="molecule-image-container">
                <img :src="predictionResult.soluteImage" alt="溶质结构" class="molecule-image" />
                <div class="image-label">溶质结构</div>
              </div>
            </div>
          </div>
        </div>
      </el-card>
    </div>
  </div>
</template>

<script setup>
import { ref, onMounted } from 'vue'
import { ElMessage } from 'element-plus'
import axios from 'axios'

const isLoggedIn = ref(false)
const predictionForm = ref({
  solventSmiles: '',
  soluteSmiles: ''
})
const predictionResult = ref(null)

// 检查登录状态
const checkLoginStatus = () => {
  const token = localStorage.getItem('token')
  isLoggedIn.value = !!token
  if (token) {
    axios.defaults.headers.common['Authorization'] = `Bearer ${token}`
  }
}

// 处理预测
const handlePredict = async () => {
  if (!predictionForm.value.solventSmiles || !predictionForm.value.soluteSmiles) {
    ElMessage.warning('请输入完整的SMILES分子式')
    return
  }

  try {
    // TODO: 这里需要添加实际的API调用
    // 临时模拟数据
    predictionResult.value = {
      score: 0.85,
      explanation: '根据分子结构分析，该溶剂和溶质具有较好的相容性。主要原因是：\n1. 分子极性相似\n2. 分子大小匹配\n3. 官能团相互作用良好',
      solventImage: 'https://via.placeholder.com/200x200?text=Solvent',
      soluteImage: 'https://via.placeholder.com/200x200?text=Solute'
    }
  } catch (error) {
    console.error('Prediction error:', error)
    ElMessage.error('预测失败：' + (error.response?.data?.detail || error.message))
  }
}

onMounted(() => {
  checkLoginStatus()
})
</script>

<style scoped>
.prediction-container {
  padding: 20px;
  max-width: 1200px;
  margin: 0 auto;
}

.input-section {
  background: white;
  padding: 20px;
  border-radius: 8px;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
  margin-bottom: 20px;
}

.result-section {
  margin-top: 20px;
}

.result-card {
  margin-bottom: 20px;
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.result-content {
  padding: 20px;
}

.result-item {
  margin-bottom: 20px;
}

.label {
  font-weight: bold;
  color: #2c785c;
  margin-right: 10px;
}

.value {
  color: #333;
}

.molecule-images {
  display: flex;
  justify-content: space-around;
  margin-top: 20px;
}

.molecule-image-container {
  text-align: center;
}

.molecule-image {
  max-width: 200px;
  max-height: 200px;
  border-radius: 8px;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

.image-label {
  margin-top: 10px;
  color: #666;
}
</style> 