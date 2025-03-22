<template>
  <div class="home-container">
    <!-- 导航栏 -->
    <el-header class="header">
      <div class="nav-container">
        <div class="logo">
          分子数据库系统
        </div>
        <div class="search-box">
          <el-input
            v-model="searchQuery"
            placeholder="请输入SMILES分子式"
            class="search-input"
            :disabled="!isLoggedIn"
          >
            <template #append>
              <el-button @click="handleSearch" :disabled="!isLoggedIn">
                <el-icon><Search /></el-icon>
              </el-button>
            </template>
          </el-input>
        </div>
        <div class="auth-buttons">
          <template v-if="!isLoggedIn">
            <el-button type="success" class="custom-button" @click="$router.push('/login')">登录</el-button>
            <el-button type="success" class="custom-button" @click="$router.push('/register')">注册</el-button>
          </template>
          <template v-else>
            <el-button type="success" class="custom-button" @click="handleLogout">退出登录</el-button>
          </template>
        </div>
      </div>
    </el-header>

    <!-- 主要内容区域 -->
    <el-main>
      <el-row :gutter="20">
        <el-col :span="18">
          <el-table
            v-loading="loading"
            :data="smilesList"
            style="width: 100%"
            border
            stripe
            highlight-current-row
          >
            <el-table-column 
              type="index" 
              label="序号" 
              width="80" 
              align="center"
            />
            <el-table-column 
              prop="smiles" 
              label="SMILES分子式" 
              min-width="300"
            >
              <template #default="scope">
                <div 
                  class="clickable-cell"
                  @click="handleRowClick(scope.row)"
                >
                  {{ scope.row.smiles }}
                </div>
              </template>
            </el-table-column>
          </el-table>
        </el-col>
        
        <el-col :span="6" v-if="isLoggedIn">
          <el-card class="history-card">
            <template #header>
              <div class="card-header">
                <span>查询历史</span>
                <el-button type="text" @click="fetchQueryHistory">刷新</el-button>
              </div>
            </template>
            <el-timeline>
              <el-timeline-item
                v-for="item in queryHistory"
                :key="item.created_at"
                :timestamp="item.created_at"
                placement="top"
              >
                <div class="history-item" @click="searchQuery = item.smiles; handleSearch()">
                  {{ item.smiles }}
                </div>
              </el-timeline-item>
            </el-timeline>
          </el-card>
        </el-col>
      </el-row>
      
      <!-- 分子图片对话框 -->
      <el-dialog
        v-model="dialogVisible"
        title="分子结构图"
        width="50%"
        center
      >
        <div class="image-container">
          <img :src="currentMoleculeImage" alt="分子结构图" v-if="currentMoleculeImage" />
        </div>
      </el-dialog>
    </el-main>
  </div>
</template>

<script setup>
import { ref, onMounted, watch } from 'vue'
import { useRouter } from 'vue-router'
import { ElMessage } from 'element-plus'
import { Search } from '@element-plus/icons-vue'
import axios from 'axios'

const router = useRouter()
const loading = ref(false)
const searchQuery = ref('')
const smilesList = ref([])
const dialogVisible = ref(false)
const currentMoleculeImage = ref('')
const isLoggedIn = ref(false)
const queryHistory = ref([])

// 检查登录状态
const checkLoginStatus = () => {
  const token = localStorage.getItem('token')
  isLoggedIn.value = !!token
  if (token) {
    // 验证 token 是否有效
    axios.get('http://localhost:8000/api/query_history/', {
      headers: {
        'Authorization': `Bearer ${token}`
      }
    }).catch(error => {
      if (error.response?.status === 401) {
        localStorage.removeItem('token')
        isLoggedIn.value = false
        router.push('/login')
      }
    })
  }
}

// 获取SMILES列表
const fetchSmilesList = async () => {
  try {
    loading.value = true
    const response = await axios.get('http://localhost:8000/api/smiles/')
    console.log('API Response:', response.data) // 添加调试日志
    
    // 处理后端返回的数据
    if (response.data && response.data.prefixes) {
      // 将 prefixes 数组转换为表格需要的格式
      smilesList.value = response.data.prefixes.map((smiles, index) => ({
        id: index + 1,
        smiles: smiles,
        created_at: new Date().toISOString() // 如果后端没有提供时间，使用当前时间
      }))
    } else {
      smilesList.value = []
      ElMessage.warning('数据格式不正确')
    }
  } catch (error) {
    console.error('Error fetching data:', error) // 添加错误日志
    ElMessage.error('获取数据失败：' + (error.response?.data?.detail || error.message))
    smilesList.value = [] // 确保在错误时设置为空数组
  } finally {
    loading.value = false
  }
}

// 处理搜索
const handleSearch = async () => {
  if (!isLoggedIn.value) {
    ElMessage.warning('请先登录后再进行搜索')
    return
  }
  
  if (!searchQuery.value) {
    ElMessage.warning('请输入要搜索的SMILES分子式')
    return
  }

  try {
    loading.value = true
    
    const token = localStorage.getItem('token')
    if (!token) {
      ElMessage.warning('登录已过期，请重新登录')
      router.push('/login')
      return
    }

    // 修改请求头和响应类型
    const response = await axios.get(`http://localhost:8000/api/query_smiles/?smiles=${encodeURIComponent(searchQuery.value)}`, {
      headers: {
        'Authorization': `Bearer ${token}`,
        'Accept': '*/*'  // 修改为接受任何类型的响应
      },
      responseType: 'blob'  // 保持 blob 类型
    })
    
    if (response.data) {
      // 根据实际的响应类型创建 blob
      const contentType = response.headers['content-type'] || 'image/png'
      const blob = new Blob([response.data], { type: contentType })
      const imageUrl = URL.createObjectURL(blob)
      currentMoleculeImage.value = imageUrl
      dialogVisible.value = true
      
      // 清理资源
      const unwatch = watch(dialogVisible, (newVal) => {
        if (!newVal) {
          URL.revokeObjectURL(imageUrl)
          unwatch()
        }
      })

      // 搜索完成后刷新历史记录
      await fetchQueryHistory()
    }
  } catch (error) {
    console.error('Search error:', error)
    if (error.response?.status === 401) {
      localStorage.removeItem('token')
      isLoggedIn.value = false
      ElMessage.error('登录已过期，请重新登录')
      router.push('/login')
    } else {
      ElMessage.error('搜索失败：' + (error.response?.data?.detail || error.message))
    }
  } finally {
    loading.value = false
  }
}

// 处理行点击
const handleRowClick = async (row) => {
  console.log('Row clicked:', row)
  
  if (!isLoggedIn.value) {
    ElMessage.warning('请先登录后再查看分子结构图')
    return
  }
  
  try {
    loading.value = true
    const encodedSmiles = encodeURIComponent(row.smiles)
    
    const response = await axios.get(`http://localhost:8000/api/structure/?smiles=${encodedSmiles}`, {
      headers: {
        'Accept': '*/*'  // 修改为接受任何类型的响应
      },
      responseType: 'blob'
    })
    
    const contentType = response.headers['content-type'] || 'image/png'
    const blob = new Blob([response.data], { type: contentType })
    const imageUrl = URL.createObjectURL(blob)
    
    currentMoleculeImage.value = imageUrl
    dialogVisible.value = true
    
    const unwatch = watch(dialogVisible, (newVal) => {
      if (!newVal) {
        URL.revokeObjectURL(imageUrl)
        unwatch()
      }
    })
    
  } catch (error) {
    console.error('Error fetching image:', error)
    ElMessage.error('获取分子结构图失败')
  } finally {
    loading.value = false
  }
}

// 处理退出登录
const handleLogout = () => {
  localStorage.removeItem('token')
  delete axios.defaults.headers.common['Authorization']
  isLoggedIn.value = false
  ElMessage.success('已退出登录')
  router.push('/login')
}

// 获取查询历史
const fetchQueryHistory = async () => {
  if (!isLoggedIn.value) return
  
  try {
    const response = await axios.get('http://localhost:8000/api/query_history/')
    queryHistory.value = response.data.history
  } catch (error) {
    console.error('Error fetching history:', error)
    ElMessage.error('获取查询历史失败')
  }
}

onMounted(() => {
  checkLoginStatus()
  fetchSmilesList()
  fetchQueryHistory()
})
</script>

<style scoped>
.home-container {
  min-height: 100vh;
  background-color: #f5f5f5;
}

.header {
  background: #2c785c;
  padding: 0;
  box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
}

.nav-container {
  height: 60px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0 20px;
}

.logo {
  color: white;
  font-size: 20px;
  font-weight: bold;
}

.search-box {
  width: 400px;
}

.search-input :deep(.el-input__wrapper) {
  background-color: rgba(255, 255, 255, 0.9);
}

.auth-buttons {
  display: flex;
  gap: 10px;
}

.custom-button {
  background-color: #84fab0 !important;
  border-color: #84fab0 !important;
  color: #2c785c !important;
}

.custom-button:hover {
  background-color: white !important;
  border-color: white !important;
  color: #2c785c !important;
}

.el-main {
  padding: 20px;
  background-color: #f5f5f5;
}

.el-table {
  margin-top: 20px;
  border-radius: 8px;
  overflow: hidden;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

:deep(.el-table__row) {
  cursor: pointer;
}

:deep(.el-table__row:hover) {
  background-color: #f0f9eb !important;
}

.image-container {
  display: flex;
  justify-content: center;
  align-items: center;
  padding: 20px;
}

.image-container img {
  max-width: 100%;
  max-height: 500px;
  border-radius: 8px;
  box-shadow: 0 2px 12px 0 rgba(0, 0, 0, 0.1);
}

.clickable-cell {
  cursor: pointer;
  padding: 5px 0;
}

.clickable-cell:hover {
  color: #2c785c;
}

.history-card {
  margin-bottom: 20px;
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.history-item {
  cursor: pointer;
  padding: 4px 0;
}

.history-item:hover {
  color: #2c785c;
}

.el-timeline {
  max-height: 600px;
  overflow-y: auto;
}
</style> 